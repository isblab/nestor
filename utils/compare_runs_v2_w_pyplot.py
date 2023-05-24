import os
import sys
import math
import yaml
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt


runs_to_compare = sys.argv[1:]
transparency = 0.75


# Start with setting default parameters for matplotlib
mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.size"] = 12


def get_all_results(all_runs: list) -> dict:
    # Read all nestor output files
    results = {}
    for run in all_runs:
        with open(os.path.join(run, "nestor_output.yaml"), "r") as resf:
            result = yaml.safe_load(resf)
        results[run.split("/")[-1]] = result
    return results


def mean_type_plotter(results: dict, key: str, ylabel: str):
    # Plots the mean values for the key
    data = []
    for parent in results:  # parent is a trial set (init_x or x_fpi)
        x_vals = []
        y_vals = []
        for run_set in results[parent]:  # runset is res_01
            all_vals = []
            for run in results[parent][run_set]:
                try:
                    val = float(results[parent][run_set][run][key])
                except ValueError as err:
                    print(f"Terminating due to the following error...\n{err}")
                    return None
                all_vals.append(val)
            x_vals.append(run_set)
            y_vals.append(np.mean(all_vals))
        data.append((x_vals, y_vals, parent))

    fig = plt.figure()
    for datum in data:
        plt.scatter(datum[0], datum[1], label=datum[2], alpha=transparency)

    plt.xlabel("Representation")
    plt.ylabel(ylabel)
    fig.legend(bbox_to_anchor=(0.9, 1.05), loc="upper right")
    fig.savefig(f"{ylabel}_comparison.png", bbox_inches="tight")
    plt.close()


def errorbar_type_plotter(results: dict, key: str, ylabel: str):
    data = []
    for parent in results:  # parent is a trial set (init_x or x_fpi)
        xvals = []
        yvals = []
        yerr = []
        for run_set in results[parent]:  # runset is res_01
            xvals.append(run_set)
            all_vals = [
                float(results[parent][run_set][run][key])
                for run in results[parent][run_set]
            ]
            yvals.append(np.mean(all_vals))
            yerr.append(np.std(all_vals) / (math.sqrt(len(all_vals))))

        data.append((xvals, yvals, yerr, parent))

    fig = plt.figure()
    for datum in data:
        plt.errorbar(
            datum[0],
            datum[1],
            yerr=datum[2],
            label=datum[3],
            fmt="o",
            alpha=transparency,
        )
    plt.xlabel("Representation")
    plt.ylabel(ylabel)
    fig.legend(bbox_to_anchor=(0.9, 1.05), loc="upper right")
    fig.savefig(f"{ylabel}_comparison.png", bbox_inches="tight")
    plt.close()


def plot_sterr(results: dict):
    """Plots standard error comparison"""
    data = []
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
        data.append((x_vals, y_vals, parent))

    fig = plt.figure()
    for datum in data:
        plt.scatter(datum[0], datum[1], label=datum[2], alpha=transparency)

    plt.xlabel("Representation")
    plt.ylabel("Log evidence")
    fig.legend(bbox_to_anchor=(0.9, 1.05), loc="upper right")
    fig.savefig("stderr_comparison.png", bbox_inches="tight")
    plt.close()


# --------------------------------------------------------------------------------------------------
# ---------------------------------------------- Main ----------------------------------------------
# --------------------------------------------------------------------------------------------------

nestor_results = get_all_results(runs_to_compare)

# plot_sterr(nestor_results, figid=figid)

toPlot_meanType: dict = {
    "nestor_process_time": "Mean NestOR process time",
    "analytical_uncertainty": "Mean analytical uncertainties",
}

toPlot_errorbarType: dict = {
    "last_iter": "Mean terations",
    "log_estimated_evidence": "Mean log(Z)",
}

for key, y_lbl in toPlot_meanType.items():
    mean_type_plotter(
        nestor_results,
        key=key,
        ylabel=y_lbl,
    )

for key, y_lbl in toPlot_errorbarType.items():
    errorbar_type_plotter(
        nestor_results,
        key=key,
        ylabel=y_lbl,
    )
