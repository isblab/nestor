from __future__ import print_function, division
import os
import glob
import time
import math
import yaml
import pickle
import RMF
import IMP
import IMP.pmi
import IMP.rmf
import IMP.isd
import IMP.pmi.dof
import numpy as np
from mpi4py import MPI


class NestedSampling:
    def __init__(self, h_param_file, nestor_restraints, rex_macro, exit_code):
        with open(h_param_file, "r") as paramf:
            self.h_params = yaml.safe_load(paramf)

        self.tic = -9999.0
        self.toc = None
        self.mcmc_step_time = None
        self.num_init_frames = self.h_params["num_init_frames"]
        self.num_frames_per_iter = self.h_params["num_frames_per_iter"]
        self.nestor_niter = self.h_params["max_nestor_iter"]
        self.rex_macro = rex_macro
        self.rex_macro.nestor_restraints = nestor_restraints
        self.rex_macro.nest = True

        self.max_plateau_hits = self.h_params["max_plateau_hits"]
        self.plateau_hits = 0
        self.max_failed_iter = self.h_params["max_failed_iterations"]
        self.failed_iter = 0

        self.Xi = 1
        self.Z = 0
        self.H = 0

        self.worst_li_list = []
        self.worst_xi_list = []
        self.log_worst_li = []
        self.log_xi = []
        self.xi = []

        self.comm_obj = MPI.COMM_WORLD
        self.finished = False
        self.termination_mode = "None"
        self.exit_code = exit_code
        self.return_vals = {}

    def sample_initial_frames(self):
        self.rex_macro.vars["number_of_frames"] = self.num_init_frames
        start_time = time.time()
        self.rex_macro.execute_macro()
        end_time = time.time()
        per_frame_sampling_time = (end_time - start_time) / self.num_init_frames
        self.mcmc_step_time = (
            per_frame_sampling_time / self.rex_macro.vars["monte_carlo_steps"]
        )

    def parse_likelihoods(self, iteration, fhead="likelihoods_"):
        sampled_likelihoods = []
        all_likelihood_binaries = glob.glob(f"{fhead}*")

        for binfile in all_likelihood_binaries:
            likelihoods = []
            with open(binfile, "rb") as rlif:
                likelihoods = pickle.load(rlif)
                for li in likelihoods:
                    sampled_likelihoods.append(li)
            os.remove(binfile)

        is_nan = False
        for li in sampled_likelihoods:
            if math.isnan(li):
                is_nan = True
        if is_nan:
            self.termination_mode = "Error: Nan found"
            self.exit_code = 11
            print("NaN found. Terminating...")
            for li in sampled_likelihoods:
                print(li)
            self.terminator(
                iteration=iteration,
                plateau_hits=self.plateau_hits,
                failed_iter=self.failed_iter,
            )

        return sampled_likelihoods

    def check_plateau(self):
        """
        Check if Li/Xi is plateuing  for consecutive samples, stop
        """

        previous_Li = self.worst_li_list[-2]
        current_Li = self.worst_li_list[-1]
        previous_Xi = self.worst_xi_list[-2]
        current_Xi = self.worst_xi_list[-1]

        if (current_Li / previous_Li) < (previous_Xi / current_Xi):
            self.plateau_hits += 1
            print(
                f"{'---'*20}\nPlateau detector hits: {self.plateau_hits}/{self.max_plateau_hits}"
            )
        else:
            self.plateau_hits = 0

        if self.plateau_hits == self.max_plateau_hits:
            self.termination_mode = "MaxPlateauHits"
            self.finished = True

    def terminator(self, iteration, plateau_hits, failed_iter):
        from math import log

        self.toc = time.time()

        if not "error" in self.termination_mode.lower():
            print(f"Estimated evidence sampled: {self.Z}")
            self.exit_code = 0
            try:
                ana_unc = math.sqrt(self.H / self.num_init_frames)
            except ValueError:
                ana_unc = "Did not compute. H was negative"
                print("Math domain error")
                self.exit_code = 13

            from matplotlib import pyplot as plt

            fig, ax = plt.subplots(1)
            ax.set_xlabel("log(Xi)")
            ax.set_ylabel("log(Li)")
            ax.clear()
            ax.plot(self.log_xi, self.log_worst_li)
            fig.savefig("log_lixi.png")
            plt.close()

            fig, ax = plt.subplots(1)
            ax.plot(self.xi, self.log_worst_li)
            ax.set_xlabel("Xi")
            ax.set_ylabel("log(Li)")
            fig.savefig("lixi.png")
            plt.close()

            self.return_vals["last_iter"] = iteration
            self.return_vals["plateau_hits"] = plateau_hits
            self.return_vals["failed_iter"] = failed_iter
            self.return_vals["obtained_information"] = self.H
            self.return_vals["analytical_uncertainty"] = ana_unc
            self.return_vals["nestor_process_time"] = self.toc - self.tic
            self.return_vals["mcmc_step_time"] = self.mcmc_step_time
            self.return_vals["log_estimated_evidence"] = log(self.Z)

        else:
            self.return_vals["run_params"] = self.h_params

        self.return_vals["termination_mode"] = self.termination_mode
        self.return_vals["exit_code"] = self.exit_code

    def compute_evidence_H(self, iteration, curr_li):
        # compute Z
        curr_xi = math.exp(-iteration / self.num_init_frames)
        curr_wi = self.Xi - curr_xi
        prev_zi = self.Z
        self.Z += curr_li * curr_wi
        curr_zi = self.Z
        self.Xi = curr_xi

        # compute H
        if iteration > 1:
            first_term = ((curr_li * curr_wi) / curr_zi) * math.log(curr_li)
            second_term = (prev_zi / curr_zi) * (self.H + math.log(prev_zi))
            self.H = first_term + second_term - math.log(curr_zi)

    def execute_nested_sampling2(self):
        self.tic = time.time()
        import matplotlib.pyplot as plt

        i = 0
        true_iter = 0
        base_process = self.comm_obj.Get_rank() == 0
        self.comm_obj.Barrier()

        if "shuffle_config.err" in os.listdir("./"):
            self.exit_code = 11

        self.comm_obj.Barrier()

        print(
            f"Exit code from the macros after communication: {self.exit_code} at rank: {self.comm_obj.Get_rank()}"
        )

        if self.exit_code is None:
            # Check for nan through small test run
            self.comm_obj.Barrier()
            if base_process:
                print(
                    f"{'-'*50}\nTest run complete, no NaN found. Continuing...\n{'-'*50}\n\n"
                )
            self.comm_obj.Barrier()

            self.sample_initial_frames()
            self.comm_obj.Barrier()

            if base_process:
                self.likelihoods = self.parse_likelihoods(iteration=true_iter)
            self.comm_obj.Barrier()

            self.rex_macro.vars["number_of_frames"] = self.num_frames_per_iter
            self.rex_macro.vars["replica_exchange_swap"] = True

            while true_iter < self.nestor_niter:
                self.comm_obj.Barrier()
                self.finished = self.comm_obj.bcast(self.finished, root=0)
                self.exit_code = self.comm_obj.bcast(self.exit_code, root=0)

                if self.exit_code is not None:
                    #'run.log' in os.listdir('./') or 'error.log' in os.listdir('./'):
                    # run log will exist if
                    # a. parse_likelihoods had a nan error in the test iter, called terminator.
                    # b. convergence criterion plateau reached, called terminator.
                    # c. convergence criterion max_failed_iterations reached, called terminator.
                    break

                if not self.finished:
                    # Other processes should not sample more models as convergence criteria i.e.
                    # a. max_failed_iterations or b. plateau triggered and the likelihoods list is
                    # unraveled to accumulate Z/H.
                    self.rex_macro.execute_macro()
                self.comm_obj.Barrier()

                if base_process:
                    if len(self.likelihoods) != 0:
                        Li = min(self.likelihoods)
                        if not self.finished:
                            newly_sampled_likelihoods = self.parse_likelihoods(
                                iteration=true_iter
                            )
                            candidate_li = max(newly_sampled_likelihoods)
                        else:  # unraveling
                            candidate_li = Li

                        if candidate_li >= Li:
                            self.likelihoods.remove(Li)

                            if not self.finished:
                                self.likelihoods.append(candidate_li)

                            # print(self.likelihoods, "\n", Li)
                            self.compute_evidence_H(iteration=i, curr_li=Li)
                            self.log_worst_li.append(math.log(Li))

                            self.log_xi.append(math.log(self.Xi))
                            self.xi.append(self.Xi)
                            self.worst_li_list.append(Li)
                            self.worst_xi_list.append(self.Xi)

                            if not self.finished:
                                if i > 1:
                                    self.check_plateau()
                                self.failed_iter = 0
                            i += 1

                        else:
                            self.failed_iter += 1
                            if self.failed_iter == self.max_failed_iter:
                                self.termination_mode = "MaxFailedIterations"
                                self.finished = True

                        true_iter += 1
                        print(
                            f'\n-----> True iteration: {true_iter} {" "*5} Calculation iteration: {i} {" "*5} Failed iteration: {self.failed_iter} {" "*5} Evidence: {self.Z} {" "*5} Terminating: {self.finished}\n'
                        )
                        if true_iter % 10 == 0:
                            from math import log

                            tempout = {
                                "True iteration": true_iter,
                                "Calculation iteration": i,
                                "Failed iteration": self.failed_iter,
                                "Log Evidence": log(self.Z),
                                "Plateau hits": self.plateau_hits,
                            }
                            with open("temporary_output.yaml", "w") as tof:
                                yaml.dump(tempout, tof)

                    else:
                        self.terminator(
                            iteration=true_iter,
                            plateau_hits=self.plateau_hits,
                            failed_iter=self.failed_iter,
                        )

                    live_fig, live_ax = plt.subplots(1)
                    live_ax.set_xlabel("log(Xi)")
                    live_ax.set_ylabel("log(Li)")
                    live_ax.plot(self.log_xi, self.log_worst_li)
                    live_fig.savefig("live_loglixi.png")
                    plt.close()

                self.comm_obj.Barrier()
                true_iter = self.comm_obj.bcast(true_iter, root=0)
                if true_iter == self.nestor_niter:
                    self.termination_mode = (
                        "Error: MaxIterations reached without convergence criteria"
                    )
                    self.exit_code = 12
                    self.exit_code = self.comm_obj.bcast(self.exit_code, root=0)
                    self.terminator(
                        iteration=true_iter,
                        plateau_hits=self.plateau_hits,
                        failed_iter=self.failed_iter,
                    )

        else:
            if base_process:
                self.termination_mode = "Error: Shuffle configuration error"
                self.exit_code = 11
                self.exit_code = self.comm_obj.bcast(self.exit_code, root=0)
                self.terminator(iteration=0, plateau_hits=0, failed_iter=0)
        self.comm_obj.Barrier()

        self.exit_code = self.comm_obj.bcast(self.exit_code, root=0)

        if base_process:
            return self.return_vals, self.exit_code
        else:
            return None, None
