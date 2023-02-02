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


def mean_type_plotter(results: dict, figid: int, key: str, ylabel: str):
    #! Plots the mean values for the key
    plt.figure(figid)
    for parent in results:  # parent is a trial set (init_x or x_fpi)
        x_vals = []
        y_vals = []
        for run_set in results[parent]:  # runset is res_01
            all_vals = []
            for run in results[parent][run_set]:
                all_vals.append(float(results[parent][run_set][run][key]))
            x_vals.append(run_set)
            y_vals.append(np.mean(all_vals))
        plt.scatter(x_vals, y_vals, label=parent, color=f"C{sys.argv.index(parent)-1}")
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(f"{ylabel}_comparison.png")


def plot_sterr_log_evi(results: dict, figid: int):
    #! Plots mean log evidence with sterr comparison
    plt.figure(figid)
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


def plot_sterr(results: dict, figid: int):
    #! Plots standard error comparison
    plt.figure(figid)
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


nestor_results = get_all_results(runs_to_compare)
figid = 0

figid += 1
plot_sterr_log_evi(nestor_results, figid=figid)
figid += 1
plot_sterr(nestor_results, figid=figid)

toPlot_meanType: dict[str:str] = {
    "nestor_process_time": "Mean NestOR process time",
    "last_iter": "Mean terations",
    "analytical_uncertainty": "Mean analytical uncertainties",
}

for key, y_lbl in toPlot_meanType.items():
    figid += 1
    mean_type_plotter(
        nestor_results,
        figid=figid,
        key=key,
        ylabel=y_lbl,
    )
