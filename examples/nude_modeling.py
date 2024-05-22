##########################################################################################
######################### IMP Modeling Script for NuDe Complex ###########################
##########################################################################################


# Imports
from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import ihm.cross_linkers
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.dof
from IMP.nestor import NestedSampling
import IMP.atom
import contextlib

# import IMP.saxs
import os
import sys

IMP.setup_from_argv(sys.argv, "Application of NestOR to the NuDe subcomplex")

dat_dir = IMP.nestor.get_example_path("input")  
run_output_dir = "run_" + "0"  # sys.argv[1]
topology_file = os.path.join(dat_dir, "topology50.txt")  # sys.argv[2]
h_param_file = os.path.join(dat_dir, "nestor_params_optrep.yaml")  # sys.argv[3]

max_shuffle_core = 5
max_shuffle_set2 = 50
rex_max_temp = 2.4

# Identify data files
sampling_adh_xl_data = os.path.join(dat_dir, "xlms/sampling_filtered_adh.dat")
sampling_bs3dss_xl_data = os.path.join(dat_dir, "xlms/sampling_filtered_bs3dss.dat")
sampling_dmtmm_xl_data = os.path.join(dat_dir, "xlms/sampling_filtered_dmtmm.dat")

evi_adh_xl_data = os.path.join(dat_dir, "xlms/evicalc_filtered_adh.dat")
evi_bs3dss_xl_data = os.path.join(dat_dir, "xlms/evicalc_filtered_adh.dat")
evi_dmtmm_xl_data = os.path.join(dat_dir, "xlms/evicalc_filtered_adh.dat")

gmm_data = os.path.join(dat_dir, "gmm/emd_22904.gmm.40.txt")
# Restraint weights
intra_xl_weight = 1.0
inter_xl_weight = 10.0
xl_weight = 10
em_weight = 1000.0


def modeling(output_dir, topology_file, h_param_file):
    mdl = IMP.Model()
    t = IMP.pmi.topology.TopologyReader(
        topology_file,
        pdb_dir=os.path.join(dat_dir, 'pdb'),
        fasta_dir=os.path.join(dat_dir, 'fasta'),
        gmm_dir=os.path.join(dat_dir, 'gmm'))
    bs = IMP.pmi.macros.BuildSystem(mdl)
    bs.add_state(t)

    root_hier, dof = bs.execute_macro(
        max_rb_trans=1,
        max_rb_rot=0.1,
        max_bead_trans=3.2,
        max_srb_trans=0.01,
        max_srb_rot=0.04,
    )

    # Fixing particles
    set1_core = []
    for prot in ["MTA1"]:
        for cp in [0, 1]:
            set1_core += IMP.atom.Selection(
                root_hier, molecule=prot, copy_index=cp, residue_indexes=range(165, 334)
            ).get_selected_particles()
    for prot in ["HDAC1"]:
        for cp in [0, 1]:
            set1_core += IMP.atom.Selection(
                root_hier, molecule=prot, copy_index=cp, residue_indexes=range(8, 377)
            ).get_selected_particles()

    set2 = []
    for prot in ["MTA1"]:
        for cp in [0, 1]:
            set2 += IMP.atom.Selection(
                root_hier, molecule=prot, copy_index=cp, residue_indexes=range(1, 165)
            ).get_selected_particles()
    for prot in ["MTA1"]:
        for cp in [0, 1]:
            set2 += IMP.atom.Selection(
                root_hier, molecule=prot, copy_index=cp, residue_indexes=range(334, 432)
            ).get_selected_particles()
    for prot in ["HDAC1"]:
        for cp in [0, 1]:
            set2 += IMP.atom.Selection(
                root_hier, molecule=prot, copy_index=cp, residue_indexes=range(1, 8)
            ).get_selected_particles()
    for prot in ["HDAC1"]:
        for cp in [0, 1]:
            set2 += IMP.atom.Selection(
                root_hier, molecule=prot, copy_index=cp, residue_indexes=range(377, 483)
            ).get_selected_particles()
    for prot in ["MBD3"]:
        set2 += IMP.atom.Selection(
            root_hier, molecule=prot, copy_index=0, residue_indexes=range(1, 296)
        ).get_selected_particles()
    for prot in ["P66A"]:
        set2 += IMP.atom.Selection(
            root_hier, molecule=prot, copy_index=0, residue_indexes=range(136, 179)
        ).get_selected_particles()
    for prot in ["RBBP4"]:
        for cp in [0, 1, 2, 3]:
            set2 += IMP.atom.Selection(
                root_hier, molecule=prot, copy_index=cp, residue_indexes=range(1, 425)
            ).get_selected_particles()

    fixed_set1_core_beads, fixed_set1_core = dof.disable_movers(
        set1_core, [IMP.core.RigidBodyMover, IMP.pmi.TransformMover]
    )

    fixed_set2_beads, fixed_set2 = dof.disable_movers(
        set2, [IMP.core.RigidBodyMover, IMP.pmi.TransformMover]
    )
    molecules = t.get_components()

    #####################################################
    ##################### RESTRAINTS ####################
    #####################################################

    output_objects = []
    nestor_restraints = []

    # -----------------------------
    # CONNECTIVITY RESTRAINT

    for m in root_hier.get_children()[0].get_children():
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
        cr.add_to_model()
        output_objects.append(cr)

    print("Connectivity restraint applied")

    # -----------------------------
    # EXCLUDED VOLUME RESTRAINT

    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=[root_hier], resolution=1000
    )
    output_objects.append(evr)

    print("Excluded volume restraint applied")

    # -------------------------
    # CROSSLINKING RESTRAINT
    # Sampling XL restraints

    xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    xldbkc.set_standard_keys()

    xldb_adh_sampling = IMP.pmi.io.crosslink.CrossLinkDataBase()
    xldb_adh_sampling.create_set_from_file(
        file_name=sampling_adh_xl_data, converter=xldbkc
    )
    xlr_adh_sampling = (
        IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xldb_adh_sampling,
            length=25,
            resolution=1,
            slope=0.0001,
            label="adh",
            weight=xl_weight,
            linker=ihm.cross_linkers.edc,
        )
    )

    xldb_bs3dss_sampling = IMP.pmi.io.crosslink.CrossLinkDataBase()
    xldb_bs3dss_sampling.create_set_from_file(
        file_name=sampling_bs3dss_xl_data, converter=xldbkc
    )
    xlr_bs3dss_sampling = (
        IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xldb_bs3dss_sampling,
            length=25,
            resolution=1,
            slope=0.0001,
            label="bs3dss",
            weight=xl_weight,
            linker=ihm.cross_linkers.bs3,
        )
    )

    xldb_dmtmm_sampling = IMP.pmi.io.crosslink.CrossLinkDataBase()
    xldb_dmtmm_sampling.create_set_from_file(
        file_name=sampling_dmtmm_xl_data, converter=xldbkc
    )
    xlr_dmtmm_sampling = (
        IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xldb_dmtmm_sampling,
            length=16,
            resolution=1,
            slope=0.0001,
            label="dmtmm",
            weight=xl_weight,
            linker=ihm.cross_linkers.dsso,
        )
    )

    output_objects.append(xlr_adh_sampling)
    output_objects.append(xlr_bs3dss_sampling)
    output_objects.append(xlr_dmtmm_sampling)

    # ----------------------------------------------------
    # Nesting XL restraints

    xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    xldbkc.set_standard_keys()

    xldb_adh_evicalc = IMP.pmi.io.crosslink.CrossLinkDataBase()
    xldb_adh_evicalc.create_set_from_file(file_name=evi_adh_xl_data, converter=xldbkc)
    xlr_adh_evicalc = (
        IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xldb_adh_evicalc,
            length=25,
            resolution=1,
            slope=0.0001,
            label="adh",
            weight=0,
            linker=ihm.cross_linkers.edc,
        )
    )

    xldb_bs3dss_evicalc = IMP.pmi.io.crosslink.CrossLinkDataBase()
    xldb_bs3dss_evicalc.create_set_from_file(
        file_name=evi_bs3dss_xl_data, converter=xldbkc
    )
    xlr_bs3dss_evicalc = (
        IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xldb_bs3dss_evicalc,
            length=25,
            resolution=1,
            slope=0.0001,
            label="bs3dss",
            weight=0,
            linker=ihm.cross_linkers.bs3,
        )
    )

    xldb_dmtmm_evicalc = IMP.pmi.io.crosslink.CrossLinkDataBase()
    xldb_dmtmm_evicalc.create_set_from_file(
        file_name=evi_dmtmm_xl_data, converter=xldbkc
    )
    xlr_dmtmm_evicalc = (
        IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
            root_hier=root_hier,
            database=xldb_dmtmm_evicalc,
            length=16,
            resolution=1,
            slope=0.0001,
            label="dmtmm",
            weight=0,
            linker=ihm.cross_linkers.dsso,
        )
    )

    output_objects.append(xlr_adh_evicalc)
    output_objects.append(xlr_bs3dss_evicalc)
    output_objects.append(xlr_dmtmm_evicalc)

    nestor_restraints.append(xlr_adh_evicalc)
    nestor_restraints.append(xlr_bs3dss_evicalc)
    nestor_restraints.append(xlr_dmtmm_evicalc)

    print("Cross-linking restraint applied")

    # EM RESTRAINT
    densities = IMP.atom.Selection(
        root_hier, representation_type=IMP.atom.DENSITIES
    ).get_selected_particles()

    emr = IMP.pmi.restraints.em.GaussianEMRestraint(
        densities,
        target_fn=gmm_data,
        slope=0.00000001,
        scale_target_to_mass=True,
        weight=0,
    )

    output_objects.append(emr)
    nestor_restraints.append(emr)
    print("EM Restraint Applied")

    #####################################################
    ###################### SAMPLING #####################
    #####################################################

    exit_code = None
    try:
        IMP.pmi.tools.shuffle_configuration(
            root_hier,
            max_translation=max_shuffle_set2,
            excluded_rigid_bodies=fixed_set1_core,
        )

        IMP.pmi.tools.shuffle_configuration(
            root_hier,
            max_translation=max_shuffle_core,
            excluded_rigid_bodies=fixed_set2,
        )

    except ValueError:
        exit_code = 11
        with open("shuffle_config.err", "w") as shufferr:
            shufferr.write(str(11))
            shufferr.flush()

    print(f"{'-'*50}\nExit code: {exit_code}\n{'-'*50}")
    dof.optimize_flexible_beads(50)

    IMP.pmi.dof.DegreesOfFreedom.enable_all_movers(dof)

    # Now, add all of the other restraints to the scoring function to start sampling
    evr.add_to_model()
    emr.add_to_model()
    xlr_adh_sampling.add_to_model()
    xlr_bs3dss_sampling.add_to_model()
    xlr_dmtmm_sampling.add_to_model()
    xlr_adh_evicalc.add_to_model()
    xlr_bs3dss_evicalc.add_to_model()
    xlr_dmtmm_evicalc.add_to_model()

    print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

    rex = IMP.pmi.macros.ReplicaExchange(
        mdl,
        root_hier=root_hier,
        
        monte_carlo_temperature=1.0,
        replica_exchange_minimum_temperature=1.0,
        replica_exchange_maximum_temperature=rex_max_temp,
        monte_carlo_sample_objects=dof.get_movers(),
        global_output_directory=run_output_dir,
        output_objects=output_objects,
        monte_carlo_steps=10,
        number_of_best_scoring_models=0,
        number_of_frames=0,
        use_nestor=True,
    )

    ns = NestedSampling(
        rex_macro=rex,
        nestor_restraints=nestor_restraints,
        h_param_file=h_param_file,
        exit_code=exit_code,
    )

    ns_return, ns_exit_code = ns.execute_nested_sampling2()
    return ns_return, ns_exit_code


with open("errors.log", "w") as errf:
    with contextlib.redirect_stdout(None):
        with contextlib.redirect_stderr(errf):
            nested_sampling_output, nested_sampling_exitcode = modeling(
                run_output_dir, topology_file, h_param_file
            )


if nested_sampling_output is not None:
    print(f"////{nested_sampling_output}")

sys.exit(nested_sampling_exitcode)
