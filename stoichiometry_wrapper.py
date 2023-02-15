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

topology = True
if "manual" in sys.argv:
    topology = False

with open(h_param_file, "r") as paramf:
    h_params = yaml.safe_load(paramf)  # type: ignore

max_allowed_runs = h_params["max_usable_threads"] // h_params["num_cores"]
parent_path = h_params["parent_dir"]

target_runs = str(h_params["num_runs"])
if "-" not in target_runs:
    target_runs = range(0, int(target_runs))
else:
    target_runs = range(int(target_runs.split("-")[0]), int(target_runs.split("-")[1]))

imp_path = " "
if "imp_path" in h_params.keys():
    imp_path = h_params["imp_path"]


###################################################
#################### Functions ####################
###################################################


def get_all_toruns(h_params) -> list[str]:
    parent_dir = h_params["parent_dir"]
    runs = []
    for stoichiometry in h_params["stoichiometries"]:
        for run in target_runs:
            run_deets = (os.path.join(parent_dir, f"{stoichiometry}"), str(run))
            runs.append(run_deets)
    return runs


def get_output(process, results: dict) -> dict:
    run_deets, proc = process
    out, _ = proc.communicate()
    result = literal_eval(out[4:])
    if run_deets[0].split("/")[-1] not in results.keys():
        results[f"{run_deets[0].split('/')[-1]}"] = {f"run_{run_deets[1]}": result}
    else:
        results[f"{run_deets[0].split('/')[-1]}"][f"run_{run_deets[1]}"] = result

    return results


def communicate_finished_proc_and_get_remaining_procs(processes, results):
    faulty_runs = []
    successful_runs = []
    terminated_runs = []

    for run_deets, proc in processes.items():
        if proc.poll() is not None:
            proc.wait()
            terminated_runs.append(run_deets)

        if proc.returncode == 11:
            faulty_runs.append(run_deets)
            shutil.rmtree(os.path.join(run_deets[0], f"run_{run_deets[1]}"))
        elif proc.returncode == 0 or proc.returncode == 13:
            successful_runs.append((run_deets, proc))

    for run in terminated_runs:
        print(
            f"Terminated: {run[0].split('/')[-1]}, run_{run[1]} with exit code: {processes[run].returncode}"
        )
        processes.pop(run)

    for p in successful_runs:
        results = get_output(p, results)

    if "nestor_output.yaml" in os.listdir(parent_path):
        with open(os.path.join(parent_path, "nestor_output.yaml"), "r") as inf:
            old_results = yaml.safe_load(inf)
            merge(results, old_results)

    with open(f"{parent_path}/nestor_output.yaml", "w") as outf:
        yaml.dump(results, outf)
        outf.flush()

    return processes, faulty_runs, successful_runs, results


def plotter(results: dict):
    all_log_z = {}
    mean_proc_time = []
    stoichiometries = []

    plt.figure(1)
    for stoichiometry in results:

        if not stoichiometry in stoichiometries:
            stoichiometries.append(int(stoichiometry.split("_")[-1]))

        log_z = []
        proc_time = []
        for _, run in results[stoichiometry].items():
            log_z.append(run["log_estimated_evidence"])
            proc_time.append(run["nestor_process_time"])

        all_log_z[stoichiometry] = log_z
        mean_proc_time.append(np.mean(proc_time))

        avg_logz = np.mean(log_z)
        stderr_logz = np.std(log_z) / math.sqrt(len(log_z))
        plt.errorbar(
            int(stoichiometry.split("_")[-1]), avg_logz, yerr=stderr_logz, fmt="o"
        )
    print(stoichiometries)
    plt.xlabel("stoichiometries")
    plt.ylabel("log(Evidence)")
    plt.savefig(
        os.path.join(
            parent_path, f"trial_{h_params['trial_name']}_evidence_errorbarplot.png"
        )
    )

    plt.figure(2)
    stoichiometries, mean_proc_time = zip(*sorted(zip(stoichiometries, mean_proc_time)))
    plt.plot(stoichiometries, mean_proc_time, c="C2", marker="o")
    plt.xlabel("Stoichiometries")
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

processes = {"Dummy": "Dummy"}
completed_runs = []
while len(list(processes.keys())) > 0:
    if "Dummy" in processes.keys():
        processes.pop("Dummy")

    if len(torun) > 0:
        curr_iter_torun = [run for run in torun]
        for stoichiometry, run_id in curr_iter_torun:
            if len(processes) < max_allowed_runs:
                if not os.path.isdir(stoichiometry):
                    os.mkdir(stoichiometry)

                os.chdir(stoichiometry)
                os.mkdir(f"run_{run_id}")
                os.chdir(f"run_{run_id}")

                if topology:
                    topf = f"topology{stoichiometry.split('/')[-1]}.txt"
                else:
                    raise Exception(
                        "NestOR for stoichiometries works with topology files only"
                    )

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
                processes[(stoichiometry, run_id)] = p
                torun.remove((stoichiometry, run_id))
                print(f"Launched: {stoichiometry.split('/')[-1]}, run_{run_id}")

            else:
                print("Waiting for free threads...")

    waiting = True
    while waiting:
        for _, p in processes.items():
            if p.poll() is not None:
                waiting = False

    (
        processes,
        faulty_runs,
        successful_runs,
        results,
    ) = communicate_finished_proc_and_get_remaining_procs(processes, results)

    for proc in successful_runs:
        completed_runs.append(proc)
    if len(processes) == 0:
        break

    if len(faulty_runs) != 0:
        for fr in faulty_runs:
            print(f"Will relaunch ({fr[0].split('/')[-1]}, run_{fr[1]})")
            torun.append(fr)


###################################################
############## Preparing the output ###############
###################################################


# for proc in completed_runs:
#     run_deets, p = proc
#     out, _ = p.communicate()
#     result = literal_eval(out[4:])

#     if run_deets[0].split("/")[-1] not in results.keys():
#         results[f"{run_deets[0].split('/')[-1]}"] = {f"run_{run_deets[1]}": result}
#     else:
#         results[f"{run_deets[0].split('/')[-1]}"][f"run_{run_deets[1]}"] = result

# if "nestor_output.yaml" in os.listdir(parent_path):
#     with open(os.path.join(parent_path, "nestor_output.yaml"), "r") as inf:
#         old_results = yaml.safe_load(inf)
#         merge(results, old_results)

# with open(f"{parent_path}/nestor_output.yaml", "w") as outf:
#     yaml.dump(results, outf)

plotter(results)
print("Done...!\n\n")
