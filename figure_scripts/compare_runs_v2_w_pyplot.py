import os
import sys
import math
import yaml
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt


TITLE = sys.argv[1]
runs_to_compare = sys.argv[2:]
transparency = 0.75
plt.rcParams["font.family"] = "sans-serif"

# Start with setting default parameters for matplotlib
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
        datum = list(datum)
        datum[0] = [str(x.split("_")[-1]) for x in datum[0]]
        plt.scatter(datum[0], datum[1], label=datum[2], alpha=transparency)

    plt.xlabel("Representation (number of residues per bead)")
    plt.ylabel(ylabel)
    plt.title(TITLE)
    # fig.legend(bbox_to_anchor=(1.15, 1.0), loc="upper right")
    fig.savefig(f"{ylabel}_comparison.png", bbox_inches="tight", dpi=600)
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
    for idx, datum in enumerate(data):
        datum = list(datum)
        datum[0] = [str(x.split("_")[-1]) for x in datum[0]]
        plt.errorbar(
            datum[0],
            datum[1],
            yerr=datum[2],
            label=datum[3],
            fmt="o",
            alpha=transparency,
            c=f"C{idx}",
        )
    plt.xlabel("Representation (number of residues per bead)")
    if "log" in ylabel:
        # plt.rcParams["text.usetex"] = True
        ylabel = "Mean log$Z$"

    plt.ylabel(ylabel)
    plt.title(TITLE)
    fig.legend(bbox_to_anchor=(1.15, 1.0), loc="upper right")
    fig.savefig(f"{ylabel}_comparison.png", bbox_inches="tight", dpi=600)
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
        plt.scatter(
            [x.split("_")[-1] for x in datum[0]],  # datum[0].split("_")[-1],
            datum[1],
            label=datum[2],
            alpha=transparency,
        )

    plt.xlabel("Representation (number of residues per bead)")
    plt.ylabel("Standard error on log(Evidence)")
    plt.title(TITLE)
    # fig.legend()  # bbox_to_anchor=(1.15, 1.0),loc="upper right"
    fig.savefig("stderr_comparison.png", bbox_inches="tight", dpi=600)
    plt.close()


# --------------------------------------------------------------------------------------------------
# ---------------------------------------------- Main ----------------------------------------------
# --------------------------------------------------------------------------------------------------

nestor_results = get_all_results(runs_to_compare)


toPlot_meanType: dict = {
    "analytical_uncertainty": "Mean analytical uncertainties",
}

toPlot_errorbarType: dict = {
    "last_iter": "Mean iterations",
    "log_estimated_evidence": "Mean log(Z)",
    "nestor_process_time": "Mean NestOR process time",
    # "mcmc_step_time": "Mean time per MCMC step",
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

plot_sterr(nestor_results)
