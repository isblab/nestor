import os
import sys
import yaml
import math
import shutil
import subprocess
import numpy as np
from mergedeep import merge
from ast import literal_eval
from matplotlib import pyplot as plt

###################################################
###################### Init #######################
###################################################

h_param_file = sys.argv[1]
topology = False

if "topology" in sys.argv:
    topology = True

with open(h_param_file, "r") as paramf:
    h_params = yaml.safe_load(paramf)  # type: ignore

max_allowed_runs = h_params["max_usable_threads"] // h_params["num_cores"]
parent_path = h_params["parent_dir"]

target_runs = str(h_params["num_runs"])
if "-" not in target_runs:
    target_runs = range(0, int(target_runs))
else:
    target_runs = range(
        int(target_runs.split("-")[0]), int(target_runs.split("-")[1]) + 1
    )

imp_path = " "
if "imp_path" in h_params.keys():
    imp_path = h_params["imp_path"]


###################################################
#################### Functions ####################
###################################################


def get_all_toruns(h_params):
    parent_dir = h_params["parent_dir"]
    runs = []
    for res in h_params["resolutions"]:
        for run in target_runs:
            run_deets = (os.path.join(parent_dir, f"res_{res}"), str(run))
            runs.append(run_deets)
    return runs


def communicate_finished_proc_and_get_remaining_procs(processes, lastcall=False):
    faulty_runs = []
    successful_runs = []
    terminated_runs = []

    for run_deets, proc in processes.items():
        if not lastcall:
            if proc.poll() is not None:
                proc.wait()
                if run_deets not in terminated_runs:
                    terminated_runs.append(run_deets)
        else:
            proc.wait()
            terminated_runs.append(run_deets)

        if proc.returncode == 11:
            faulty_runs.append(run_deets)
            shutil.rmtree(os.path.join(run_deets[0], f"run_{run_deets[1]}"))
        elif proc.returncode == 0:
            successful_runs.append((run_deets, proc))

    for run in terminated_runs:
        print(f"Terminated: {run[0].split('/')[-1]}, run_{run[1]}")
        processes.pop(run)

    return processes, faulty_runs, successful_runs


def plotter(results: dict):
    all_log_z = {}
    mean_proc_time = []
    resolutions = []

    plt.figure(1)
    for resolution in results:
        if not resolution[4:] in resolutions:
            resolutions.append(int(resolution[4:]))

        log_z = []
        proc_time = []
        for _, run in results[resolution].items():
            log_z.append(run["log_estimated_evidence"])
            proc_time.append(run["nestor_process_time"])

        all_log_z[resolution] = log_z
        mean_proc_time.append(np.mean(proc_time))

        avg_logz = np.mean(log_z)
        stderr_logz = np.std(log_z) / math.sqrt(len(log_z))
        plt.errorbar(int(resolution[4:]), avg_logz, yerr=stderr_logz, fmt="o")

    plt.xlabel("Resolutions")
    plt.ylabel("log(Evidence)")
    plt.savefig(
        os.path.join(
            parent_path, f"trial_{h_params['trial_name']}_evidence_errorbarplot.png"
        )
    )

    plt.figure(2)
    resolutions, mean_proc_time = zip(*sorted(zip(resolutions, mean_proc_time)))
    plt.plot(resolutions, mean_proc_time, c="C2", marker="o")
    plt.xlabel("Resolutions")
    plt.ylabel("Nested sampling process time")
    plt.savefig(
        os.path.join(parent_path, f"trial_{h_params['trial_name']}_proctime.png")
    )


###################################################
###################### Main #######################
###################################################

if "skip_calc" in sys.argv:
    with open(os.path.join(parent_path, "nestor_output.yaml"), "r") as outf:
        results = yaml.safe_load(outf)
        plotter(results)
        exit()

torun = get_all_toruns(h_params)
results = {}

processes = {}
completed_runs = []
while len(torun) > 0:
    curr_iter_torun = [run for run in torun]
    for res, run_id in curr_iter_torun:
        if len(processes) < max_allowed_runs:
            if not os.path.isdir(res):
                os.mkdir(res)

            os.chdir(res)
            os.mkdir(f"run_{run_id}")
            os.chdir(f"run_{run_id}")

            if topology:
                topf = f"topology{res.split('/')[-1].split('_')[-1]}.txt"
            else:
                topf = res.split("/")[-1].split("_")[-1]

            run_cmd = [
                "mpirun",
                "-n",
                str(h_params["num_cores"]),
                h_params["imp_path"],
                "python",
                h_params["modeling_script_path"],
                str(run_id),
                topf,
                h_param_file,
            ]

            p = subprocess.Popen(
                run_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            processes[(res, run_id)] = p
            torun.remove((res, run_id))
            print(f"Launched: {res.split('/')[-1]}, run_{run_id}")

        else:
            print("\nWaiting for free threads...\n")

    waiting = True
    while waiting:
        for _, p in processes.items():
            if p.poll() is not None:
                waiting = False

    (
        processes,
        faulty_runs,
        successful_runs,
    ) = communicate_finished_proc_and_get_remaining_procs(processes)

    for proc in successful_runs:
        completed_runs.append(proc)
    if len(processes) == 0:
        break

    if len(faulty_runs) != 0:
        for fr in faulty_runs:
            print(f"Will relaunch ({fr[0].split('/')[-1]}, run_{fr[1]})")
            torun.append(fr)


print(f"Waiting for {len(processes.keys())} processes to terminate...")
(
    processes,
    faulty_runs,
    successful_runs,
) = communicate_finished_proc_and_get_remaining_procs(processes, lastcall=True)


###################################################
############## Preparing the output ###############
###################################################

print("Performing housekeeping tasks")
for proc in successful_runs:
    completed_runs.append(proc)

for proc in completed_runs:
    run_deets, p = proc
    out, _ = p.communicate()

    result = literal_eval(out[4:])
    print(run_deets, result)

    if run_deets[0].split("/")[-1] not in results.keys():
        results[f"{run_deets[0].split('/')[-1]}"] = {f"run_{run_deets[1]}": result}
    else:
        results[f"{run_deets[0].split('/')[-1]}"][f"run_{run_deets[1]}"] = result

if "nestor_output.yaml" in os.listdir(parent_path):
    with open(os.path.join(parent_path, "nestor_output.yaml"), "r") as inf:
        old_results = yaml.safe_load(inf)
        merge(results, old_results)

with open(f"{parent_path}/nestor_output.yaml", "w") as outf:
    yaml.dump(results, outf)

plotter(results)
