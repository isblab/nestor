import os
import sys
import yaml
import math
import numpy as np
import matplotlib.pyplot as plt


parent_path = "/home/shreyas/Projects/cgopt/systems/to_ms/latest/optrep_nestor/"
systems = ["gtusc", "rnapolii", "mhm", "nude"]
nestor_results_paths = []
nestor_params_paths = []


for i in systems:
    nestor_results_paths.append(
        os.path.join(parent_path, i, "5runs", "nestor_output.yaml")
    )
    nestor_params_paths.append(os.path.join(parent_path, i, "nestor_params.yaml"))

for result, params in zip(nestor_results_paths, nestor_params_paths):
    with open(result, "r") as resf:
        with open(params, "r") as paramf:
            nestor_out = yaml.safe_load(resf)
            nestor_params = yaml.safe_load(paramf)

    representations = list(nestor_out.keys())
    y_vals = []
    y_sterr = []
    for rep in representations:
        times = []
        for run in list(nestor_out[rep].keys()):
            total_time = nestor_out[rep][run]["nestor_process_time"]

            num_iter = int(
                nestor_out[rep][run]["last_iter"]
            )  # ? This is the true iteration count

            fpi = int(nestor_params["num_frames_per_iter"])
            ninit_frames = int(nestor_params["num_init_frames"])

            total_frames = (fpi * num_iter) + ninit_frames
            time_per_frame = total_time / total_frames

            times.append(time_per_frame)
        y_vals.append(np.mean(times))
        y_sterr.append(np.std(times) / math.sqrt(len(times)))

    plt.errorbar(representations, y_vals, color="C2", yerr=y_sterr, marker="o")
    plt.xlabel("Representation")
    plt.ylabel("Average time required per frame (s)")

    outpath = "/".join(result.split("/")[:-1])
    plt.savefig(os.path.join(outpath, "Nestor time per frame.png"))
