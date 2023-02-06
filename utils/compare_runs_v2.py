import os
import sys
import yaml
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objects as go

runs_to_compare = sys.argv[1:]


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
    data = []
    for parent in results:  # parent is a trial set (init_x or x_fpi)
        x_vals = []
        y_vals = []
        for run_set in results[parent]:  # runset is res_01
            all_vals = []
            for run in results[parent][run_set]:
                try:
                    val = float(results[parent][run_set][run][key])
                except ValueError:
                    return None
                all_vals.append(val)
            x_vals.append(run_set)
            y_vals.append(np.mean(all_vals))
        data.append((x_vals, y_vals, parent))

    fig = go.Figure()
    for datum in data:
        fig.add_scatter(x=datum[0], y=datum[1], name=datum[2], mode="markers")

    fig.update_layout(
        xaxis_title="Resolutions",
        yaxis_title=ylabel,
        legend_title="Legend",
        font=dict(family="Arial", size=18),
        template="simple_white",
    )
    fig.write_html(f"{ylabel}_comparison.html")


def errorbar_type_plotter(results: dict, figid: int, key: str, ylabel: str):
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

        data.append((xvals, yvals, dict(type="data", array=yerr, width=2), parent))

    fig = go.Figure()
    for datum in data:
        fig.add_scatter(
            x=datum[0], y=datum[1], error_y=datum[2], name=datum[3], mode="markers"
        )

    fig.update_layout(
        xaxis_title="Resolutions",
        yaxis_title=ylabel,
        legend_title="Legend",
        font=dict(family="Arial", size=18),
        template="simple_white",
    )
    fig.write_html(f"{ylabel}_comparison.html")


def plot_sterr(results: dict, figid: int):
    #! Plots standard error comparison
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

    fig = go.Figure()
    for datum in data:
        fig.add_scatter(x=datum[0], y=datum[1], name=datum[2], mode="markers")

    fig.update_layout(
        xaxis_title="Resolutions",
        yaxis_title="Standard error on log(Evidence)",
        legend_title="Legend",
        font=dict(family="Arial", size=18),
        template="simple_white",
    )
    fig.write_html("Standard error on log(Evidence)_comparison.html")


# --------------------------------------------------------------------------------------------------
# ---------------------------------------------- Main ----------------------------------------------
# --------------------------------------------------------------------------------------------------

nestor_results = get_all_results(runs_to_compare)

figid = 1
plot_sterr(nestor_results, figid=figid)

toPlot_meanType: dict[str:str] = {
    "nestor_process_time": "Mean NestOR process time",
    "last_iter": "Mean terations",
    "analytical_uncertainty": "Mean analytical uncertainties",
}

toPlot_errorbarType: dict[str:str] = {
    "last_iter": "Mean terations",
    "log_estimated_evidence": "Mean log(Z)",
}


for key, y_lbl in toPlot_meanType.items():
    figid += 1
    mean_type_plotter(
        nestor_results,
        figid=figid,
        key=key,
        ylabel=y_lbl,
    )

for key, y_lbl in toPlot_errorbarType.items():
    figid += 1
    errorbar_type_plotter(
        nestor_results,
        figid=figid,
        key=key,
        ylabel=y_lbl,
    )
