import os
import sys
import yaml
import math
import numpy as np
import matplotlib.pyplot as plt

xlabel = sys.argv[1]
runs_to_compare = sys.argv[2:]


def get_all_results(all_runs: list) -> dict:
    #! Read all nestor output files
    nestor_results = {}
    for run in all_runs:
        with open(os.path.join(run, "nestor_output.yaml"), "r") as resf:
            result = yaml.safe_load(resf)
        nestor_results[run] = result
    return nestor_results


def plot_nestor_process_times(results: dict):
    #! Plots mean process time comparison
    plt.figure(1)
    for parent in results:  # parent is a trial set (init_x or x_fpi)
        x_vals = []
        y_vals = []
        for run_set in results[parent]:  # runset is res_01
            total_time = []
            for run in results[parent][run_set]:
                total_time.append(
                    float(results[parent][run_set][run]["nestor_process_time"])
                )
            x_vals.append(run_set)
            y_vals.append(np.mean(total_time))
        plt.scatter(x_vals, y_vals, label=parent, color=f"C{sys.argv.index(parent)-1}")
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel("Nestor process time")
    plt.savefig(f"Process time comparison.png")


def plot_sterr_log_evi(results: dict):
    #! Plots mean log evidence with sterr comparison
    plt.figure(2)
    done_runs = []
    for parent in results:  # parent is a trial set (init_x or x_fpi)
        for run_set in results[parent]:  # runset is res_01
            log_evi = []
            for run in results[parent][run_set]:
                log_evi.append(
                    float(results[parent][run_set][run]["log_estimated_evidence"])
                )
            mean_log_evi = np.mean(log_evi)
            stderr_log_evi = np.std(log_evi) / (math.sqrt(len(log_evi)))
            if parent not in done_runs:
                plt.errorbar(
                    run_set,
                    mean_log_evi,
                    yerr=stderr_log_evi,
                    fmt="o",
                    color=f"C{sys.argv.index(parent)-1}",
                    label=parent,
                )
                done_runs.append(parent)

            else:
                plt.errorbar(
                    run_set,
                    mean_log_evi,
                    yerr=stderr_log_evi,
                    fmt="o",
                    color=f"C{sys.argv.index(parent)-1}",
                )

    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel("Mean log(evidence) with standard error")
    plt.savefig(f"Mean log evidence comparison with sterr.png")


def plot_sterr(results: dict):
    #! Plots standard error comparison
    plt.figure(3)
    for parent in results:  # parent is a trial set (init_x or x_fpi)
        x_vals = []
        y_vals = []
        for run_set in results[parent]:  # runset is res_01
            log_evi = []
            for run in results[parent][run_set]:
                log_evi.append(
                    float(results[parent][run_set][run]["log_estimated_evidence"])
                )
            stderr_log_evi = np.std(log_evi) / (math.sqrt(len(log_evi)))
            x_vals.append(run_set)
            y_vals.append(stderr_log_evi)
        plt.scatter(x_vals, y_vals, label=parent, color=f"C{sys.argv.index(parent)-1}")

    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel("Standard error on log(Evidence)")
    plt.savefig(f"Standard error comparison.png")


def plot_analytical_uncertainty(results: dict):
    #! Plots analytical uncertainty comparison
    plt.figure(4)
    for parent in results:  # parent is a trial set (init_x or x_fpi)
        x_vals = []
        y_vals = []
        for run_set in results[parent]:  # runset is res_01
            ana_unc = []
            for run in results[parent][run_set]:
                ana_unc.append(
                    float(results[parent][run_set][run]["analytical_uncertainty"])
                )
            x_vals.append(run_set)
            y_vals.append(np.mean(ana_unc))
        plt.scatter(x_vals, y_vals, label=parent, color=f"C{sys.argv.index(parent)-1}")
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel("Analytical uncertainty")
    plt.savefig(f"Analytical uncertainty comparison.png")


nestor_results = get_all_results(runs_to_compare)
plot_nestor_process_times(nestor_results)
plot_sterr_log_evi(nestor_results)
plot_sterr(nestor_results)
plot_analytical_uncertainty(nestor_results)
