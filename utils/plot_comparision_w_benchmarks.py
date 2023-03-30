import os
import sys
import yaml
import math
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


########################################################################################################
############################################ Setting it  up ############################################
########################################################################################################

parent_path = "/home/shreyas/Projects/cgopt/systems/to_ms/latest/optrep_nestor/"
benchmark_results_path = "/home/shreyas/Projects/cgopt/systems/to_ms/latest/optrep_nestor/benchmarks/benchmarks_combined.yaml"

systems = ["gtusc", "rnapolii", "mhm", "nude"]
nestor_results_paths = []
for i in systems:
    nestor_results_paths.append(
        os.path.join(parent_path, i, "5runs", "nestor_output.yaml")
    )

fields = ["GMM Cross-correlation", "Model Precision", "Total XL score"]


for nestor_output_path in tqdm(nestor_results_paths):
    system_name = nestor_output_path.split("/")[-3]

    with open(benchmark_results_path, "r") as bfile:
        benchmark_results = yaml.safe_load(bfile)[system_name]
    with open(nestor_output_path, "r") as nfile:
        nestor_results = yaml.safe_load(nfile)

    for field in fields:
        # Field check
        if not field in benchmark_results[list(benchmark_results.keys())[0]]:
            break

        b_xvals, b_yvals = [], []  # benchmark values
        for representation in benchmark_results:
            b_xvals.append(representation)
            b_yvals.append(benchmark_results[representation][field])

        n_xvals, n_yvals, n_ystrr = [], [], []  # mean log evidence from nestor
        for representation in nestor_results:
            n_xvals.append(representation)
            temp_yvals = []
            for run in nestor_results[representation]:
                temp_yvals.append(
                    nestor_results[representation][run]["log_estimated_evidence"]
                )
            n_yvals.append(-1 * np.mean(temp_yvals))
            sterr = np.std(temp_yvals) / math.sqrt((len(temp_yvals)))
            n_ystrr.append(sterr)

        fig, ax1 = plt.subplots()
        ax1.scatter(
            b_xvals, b_yvals, color="C2", marker="o", label=f"Benchmark: {field}"
        )
        ax1.set_xlabel("Resolution")
        ax1.set_ylabel("Benchmark scores")

        ax2 = ax1.twinx()
        ax2.errorbar(
            n_xvals,
            n_yvals,
            yerr=n_ystrr,
            color="C1",
            fmt="o",
            label="Nestor result",
        )
        ax2.set_ylabel("Mean -log(Evidence from Nestor)")
        fig.legend()

        output_fig_path = os.path.join(
            "/".join(nestor_output_path.split("/")[:-1]),
            f"Comparision_with_benchmark - {field}.png",
        )
        fig.savefig(output_fig_path)
