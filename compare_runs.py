import os
import sys
import yaml
import math
import numpy as np
from matplotlib import pyplot as plt


results_files = sys.argv[1:]


def get_mean_log_evidence_and_sterr(runs: dict) -> tuple[float, float]:
    log_evidences = []
    for run_results in runs:
        log_evidences.append(runs[run_results]["log_estimated_evidence"])

    stderr = np.std(log_evidences) / math.sqrt(len(log_evidences))
    return (np.mean(log_evidences), stderr)


def get_mean_process_time(runs: dict) -> float:
    proc_times = []
    for run_results in runs:
        proc_times.append(runs[run_results]["nestor_process_time"])

    return np.mean(proc_times)


for results_file in results_files:
    with open(results_file, "r") as rf:
        nestor_results = yaml.safe_load(rf)

    xvals = []
    log_evidence_deets = []
    mean_proc_times = []
    for resolution in nestor_results:
        xvals.append(int(resolution.split("_")[-1]))
        log_evidence_deets.append(
            get_mean_log_evidence_and_sterr(nestor_results[resolution])
        )
        mean_proc_times.append(get_mean_process_time(nestor_results[resolution]))

    plt.figure(1)
    plt.errorbar(
        xvals,
        [i[0] for i in log_evidence_deets],
        yerr=[i[1] for i in log_evidence_deets],
        fmt="o-",
        label=results_file.split("/")[-2],
    )

    plt.figure(2)
    plt.plot(xvals, mean_proc_times, marker="o", label=results_file.split("/")[-2])

plt.figure(1)
plt.xlabel("Resolutions")
plt.ylabel("log(Evidence)")
plt.legend()
plt.savefig("tmp1.png")

plt.figure(2)
plt.xlabel("Resolutions")
plt.ylabel("Nestor process time (sec)")
plt.legend()
plt.savefig("tmp2.png")
