import os
import shutil
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
import IMP.nestor.wrapper_v6 as wrapper_v6


class Tests(IMP.test.TestCase):
    def prepare_system(self, topology_fname="topology5.txt"):
        """Set the system up for Nested Sampling runs"""
        topology_file = self.get_input_file_name(topology_fname)
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
            file_name=self.get_input_file_name("xl_dataset_test_io_crosslink_map.txt"),
            converter=xldbkc,
        )
        xlr_test = (
            IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,
                database=xldb_test,
                length=16,
                resolution=1,
                slope=0.0001,
                label="test",
                weight=0,
                linker=ihm.cross_linkers.edc,
            )
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
            h_param_file=self.get_input_file_name("nestor_params_optrep.yaml"),
            exit_code=exit_code,
        )
        return ns

    def run_base_run(self):
        """Was used to get the preliminary output for test_reproducibility() function. It is not used anymore."""
        results = {}
        for i in range(10):
            ns = self.prepare_system(topology_fname="topology5.txt")
            ns_output, _ = ns.execute_nested_sampling2()
            results[i] = ns_output
        with open(
            "input/nestor_output.yaml",
            "w",
        ) as yaml_out:
            yaml.dump(results, yaml_out)

    def check_individual_run_output_file_creation(self):
        """
        Check if all output files are created from the Nested Sampling run"""
        expected_files = [
            "live_loglixi.png",
            "lixi.png",
            "log_lixi.png",
            "nested_0.rmf3",
            "temporary_output.yaml",
        ]

        all_files_in_dir = os.listdir(os.getcwd())
        for exp_file in expected_files:
            self.assertTrue(exp_file in all_files_in_dir)

    def test_ns_initial_sampling(self):
        """Test initial sampling and likelihood parsing"""
        ns = self.prepare_system("topology5.txt")
        ns.sample_initial_frames()
        likelihoods = ns.parse_likelihoods(iteration=0)

        with open(self.get_input_file_name("nestor_params_optrep.yaml"), "r") as yamf:
            expectations = yaml.safe_load(yamf)

        self.assertEqual(len(likelihoods), int(expectations["num_init_frames"]))

    def test_reproducibility(self):
        """Check if the results are reproducible"""
        with open(self.get_input_file_name("nestor_output.yaml"), "r") as ori_out:
            expected_result = yaml.safe_load(ori_out)
            expected_logz = [
                expected_result[i]["log_estimated_evidence"] for i in expected_result
            ]

            new_results = []
            for _ in range(3):
                ns = self.prepare_system("topology5.txt")
                ns_output, _ = ns.execute_nested_sampling2()
                new_results.append(ns_output["log_estimated_evidence"])
                self.check_individual_run_output_file_creation()

            mean_expected, std_expected = np.mean(np.array(expected_logz)), np.std(
                np.array(expected_logz)
            )

            mean_new = np.mean(np.array(new_results))

            lower_bound = mean_expected - (1.96 * std_expected)
            upper_bound = mean_expected + (1.96 * std_expected)
            self.assertTrue(lower_bound <= mean_new <= upper_bound)


if __name__ == "__main__":
    IMP.test.main()
