import os, sys
import numpy as np
import yaml
import IMP
import ihm.cross_linkers
import IMP.test
from IMP.nestor import nestor
import IMP.pmi.topology
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.tools
import IMP.pmi.macros
import contextlib


def prepare_system():
    """Set the system up for Nested Sampling runs"""
    topology_file = "../../../input/topology.txt"
    mdl = IMP.Model()
    t = IMP.pmi.topology.TopologyReader(topology_file)
    bs = IMP.pmi.macros.BuildSystem(mdl)
    bs.add_state(t)

    root_hier, dof = bs.execute_macro(
        max_rb_trans=1,
        max_rb_rot=0.1,
        max_bead_trans=3.2,
        max_srb_trans=0.01,
        max_srb_rot=0.04,
    )

    nestor_restraints = []
    for m in root_hier.get_children()[0].get_children():
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
        cr.add_to_model()
    evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
        included_objects=[root_hier], resolution=1000
    )
    evr.add_to_model()

    xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    xldbkc.set_protein1_key("P1")
    xldbkc.set_protein2_key("P2")
    xldbkc.set_residue1_key("R1")
    xldbkc.set_residue2_key("R2")
    xldb_test = IMP.pmi.io.crosslink.CrossLinkDataBase()
    xldb_test.create_set_from_file(
        file_name="../../../input/xl_dataset_test_io_crosslink_map.txt",
        converter=xldbkc,
    )
    xlr_test = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=root_hier,
        database=xldb_test,
        length=16,
        resolution=1,
        slope=0.0001,
        label="test",
        weight=0,
        linker=ihm.cross_linkers.edc,
    )
    nestor_restraints.append(xlr_test)

    def shuffle(rh):
        try:
            IMP.pmi.tools.shuffle_configuration(rh, max_translation=100)
            return None
        except ValueError:
            return 11

    exit_code = 11
    while exit_code is not None:
        exit_code = shuffle(root_hier)

    xlr_test.add_to_model()

    rex = IMP.pmi.macros.ReplicaExchange(
        mdl,
        root_hier=root_hier,  # pass the root hierarchy
        monte_carlo_temperature=1.0,
        replica_exchange_minimum_temperature=1.0,
        replica_exchange_maximum_temperature=2.4,
        monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
        global_output_directory="ns_test_outputs/",  # The output directory for this sampling run.
        monte_carlo_steps=10,  # Number of MC steps between writing frames
        number_of_best_scoring_models=0,  # set >0 to store best PDB files (but this is slow)
        number_of_frames=0,  # Total number of frames to run / write to the RMF file.
        use_nestor=True,
    )
    ns = nestor.NestedSampling(
        rex_macro=rex,
        nestor_restraints=nestor_restraints,
        h_param_file="../../../input/nestor_params_optrep.yaml",
        exit_code=exit_code,
    )
    return ns


with open("errors.log", "w") as errf:
    with contextlib.redirect_stdout(None):
        with contextlib.redirect_stderr(errf):
            ns = prepare_system()
            ns_output, ns_exitcode = ns.execute_nested_sampling2()

if ns_output is not None:
    print(f"////{ns_output}")

sys.exit(ns_exitcode)
