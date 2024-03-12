import os
import sys
import math
import yaml
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

target_run: str = sys.argv[1]
output_fname: str = "sterr_evi_and_proctime"
system_name = sys.argv[2]

with open(os.path.join(target_run, "nestor_output.yaml"), "r") as target_file:
    nestor_results: dict = yaml.safe_load(target_file)

representations: list = []
log_evi_mean_sterr: list = []
proctime_mean_sterr: list = []
for k in nestor_results:
    representations.append(k.split("_")[-1])
    log_evi, proctime = [], []

    for k1 in nestor_results[k]:
        log_evi.append(nestor_results[k][k1]["log_estimated_evidence"])
        proctime.append(nestor_results[k][k1]["mcmc_step_time"])

    log_evi_mean_sterr.append(
        (np.mean(log_evi), np.std(log_evi) / math.sqrt(len(log_evi)))
    )
    proctime_mean_sterr.append(
        (np.mean(proctime), np.std(proctime) / math.sqrt(len(proctime)))
    )

log_evi_mean_sterr = np.array(log_evi_mean_sterr)
proctime_mean_sterr = np.array(proctime_mean_sterr)

fig, ax1 = plt.subplots()
ax1.errorbar(
    x=representations,
    y=log_evi_mean_sterr[:, 0],
    yerr=log_evi_mean_sterr[:, 1],
    fmt="o",
    c="dodgerblue",
    label="Log(Evidence)",
)
# plt.rcParams["text.usetex"] = True
ylabel = "Mean log$Z$"
ax1.set_ylabel(ylabel)
ax1.set_xlabel("Representations")

ax2 = plt.twinx(ax=ax1)

ax2.scatter(
    x=representations,
    y=proctime_mean_sterr[:, 0],
    # yerr=proctime_mean_sterr[:, 1],
    # fmt="o",
    c="green",
    label="NestOR process time",
)
ax2.set_ylabel("Time per MCMC sampling step (sec)")
# fig.legend()
plt.title(system_name)
fig.savefig(f"{output_fname}.png", dpi=1200)
