"""
#############################################
##  IMP Tutorial Script
##
#############################################
#
# Short modeling script combining EM and Crosslinking data
# to localize two domains of RNA Polymerase II
#
# Authors: Riccardo Pellarin, Charles Greenberg, Daniel Saltzberg
#
# References: Papers where this data is shown...
#
"""
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic

# import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.pmi.topology
import ihm.cross_linkers
import sys, os
import contextlib

# ---------------------------
# Define Input Files
# ---------------------------
datadirectory = "/home/shreyas/Projects/cgopt/systems/to_ms/ninit/rnapolii/data/"

output_dir = "run_" + sys.argv[1]
topology_file = datadirectory + sys.argv[2]
h_param_file = sys.argv[3]

target_gmm_file = datadirectory + "emd_1883.map.mrc.gmm.50.txt"


def modeling(output_dir, topology_file, h_param_file):
    m = IMP.Model()
    topology = IMP.pmi.topology.TopologyReader(
        topology_file,
        pdb_dir=datadirectory,
        fasta_dir=datadirectory,
        gmm_dir=datadirectory,
    )

    bs = IMP.pmi.macros.BuildSystem(m)
    bs.add_state(topology)
    root_hier, dof = bs.execute_macro(
        max_rb_trans=4.0,
        max_rb_rot=0.3,
        max_bead_trans=4.0,
        max_srb_trans=4.0,
        max_srb_rot=0.3,
    )

    fixed_particles = []
    for prot in [
        "Rpb1",
        "Rpb2",
        "Rpb3",
        "Rpb5",
        "Rpb6",
        "Rpb8",
        "Rpb9",
        "Rpb10",
        "Rpb11",
        "Rpb12",
    ]:
        fixed_particles += IMP.atom.Selection(
            root_hier, molecule=prot
        ).get_selected_particles()

    fixed_beads, fixed_rbs = dof.disable_movers(
        fixed_particles, [IMP.core.RigidBodyMover, IMP.pmi.TransformMover]
    )

    exit_code = None
    try:
        IMP.pmi.tools.shuffle_configuration(
            root_hier,
            excluded_rigid_bodies=fixed_rbs,
            max_translation=50,
            verbose=False,
            cutoff=5.0,
            niterations=100,
        )

    except ValueError:
        exit_code = 11

    if exit_code is None:
        # -----------------------------------
        # Define Scoring Function Components
        # -----------------------------------

        outputobjects = []  # reporter objects (for stat files)
        nester_restraints = []

        # ---- Sequence connectivity restraint
        mols = IMP.pmi.tools.get_molecules(root_hier)
        for mol in mols:
            molname = mol.get_name()
            IMP.pmi.tools.display_bonds(mol)
            cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
                mol, scale=2.0
            )
            cr.add_to_model()
            cr.set_label(molname)
            outputobjects.append(cr)

        # ---- Excluded Volume Restraint
        ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
            included_objects=root_hier, resolution=10
        )
        ev.add_to_model()
        outputobjects.append(ev)

        # ---- Electron Microscopy Restraint
        em_components = IMP.pmi.tools.get_densities(root_hier)

        gemt = IMP.pmi.restraints.em.GaussianEMRestraint(
            em_components,
            target_gmm_file,
            scale_target_to_mass=True,
            slope=0.000001,
            weight=0,
        )  # 80.0)
        gemt.add_to_model()
        outputobjects.append(gemt)
        nester_restraints.append(gemt)

        # ---- Crosslinking mass spectrometry restraint
        # Sampling XLs

        xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
        xldbkwc.set_standard_keys()

        xl1db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
        xl1db.create_set_from_file(
            f"{datadirectory}/sampling_polii_xlinks_w_header.csv"
        )

        xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xl1db,
            length=21.0,
            slope=0.01,
            resolution=1.0,
            label="Trnka",
            linker=ihm.cross_linkers.dss,
            weight=1,
        )

        xl1.add_to_model()  # crosslink must be added to the model
        outputobjects.append(xl1)

        # Nesting XLs
        xldbkwc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
        xldbkwc.set_standard_keys()

        xl2db = IMP.pmi.io.crosslink.CrossLinkDataBase(xldbkwc)
        xl2db.create_set_from_file(f"{datadirectory}/evicalc_polii_xlinks_w_header.csv")

        xl2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xl2db,
            length=21.0,
            slope=0.01,
            resolution=1.0,
            label="Chen",
            linker=ihm.cross_linkers.bs3,
            weight=0,
        )
        xl2.add_to_model()
        outputobjects.append(xl2)
        nester_restraints.append(xl2)

        for i in outputobjects:
            print(i)
        print("Nesting using these restraints:")
        for j in nester_restraints:
            print(j)
        # --------------------------
        # Monte-Carlo Sampling
        # --------------------------

        mc1 = IMP.pmi.macros.ReplicaExchange0(
            m,
            root_hier=root_hier,
            monte_carlo_sample_objects=dof.get_movers(),
            output_objects=outputobjects,
            monte_carlo_temperature=1.0,
            replica_exchange_minimum_temperature=1.0,
            replica_exchange_maximum_temperature=2.5,
            number_of_best_scoring_models=0,
            monte_carlo_steps=10,
            number_of_frames=0,
            global_output_directory=output_dir,
            use_nester=True,
        )

        ns: IMP.pmi.macros.NestedSampling = IMP.pmi.macros.NestedSampling(
            rex_macro=mc1,
            nester_restraints=nester_restraints,
            h_param_file=h_param_file,
            exit_code=None,
        )

    else:
        ns = IMP.pmi.macros.NestedSampling(
            rex_macro=None,
            nester_restraints=None,
            h_param_file=h_param_file,
            exit_code=exit_code,
        )

    ns_return, ns_exit_code = ns.execute_nested_sampling2()

    return ns_return, ns_exit_code


with contextlib.redirect_stdout(None):
    nested_sampling_output, nested_sampling_exitcode = modeling(
        output_dir, topology_file, h_param_file
    )
if nested_sampling_output is not None:
    print(f"////{nested_sampling_output}")

sys.exit(nested_sampling_exitcode)
