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
#################### Functions ####################
###################################################


def get_all_toruns(h_params, target_runs):
    parent_dir = h_params["parent_dir"]
    runs = []
    for res in h_params["resolutions"]:
        for run in target_runs:
            run_deets = (os.path.join(parent_dir, f"res_{res}"), str(run))
            runs.append(run_deets)
    return runs


def get_curr_processes_and_terminated_runs(processes: dict):
    faulty_runs = []
    successful_runs = []
    terminated_runs = []

    for run_deets, proc in processes.items():
        if proc.poll() is not None:
            proc.wait()
            terminated_runs.append(run_deets)

        if proc.returncode == 11 or proc.returncode == 12:
            faulty_runs.append((run_deets, proc))
            if proc.returncode == 11:
                shutil.rmtree(os.path.join(run_deets[0], f"run_{run_deets[1]}"))
        elif proc.returncode == 0:
            successful_runs.append((run_deets, proc))

    for run in terminated_runs:
        print(
            f"Terminated: {run[0].split('/')[-1]}, run_{run[1]} with exit code: {processes[run].returncode}"
        )
        if processes[run].returncode != 0:
            print(f"Error:\n{processes[run].stderr.read()}")

        processes.pop(run)

    return processes, faulty_runs, successful_runs


def plotter(results: dict, h_params):
    all_log_z = {}
    mean_proc_time = []
    mean_per_step_time = []
    resolutions = []

    plt.figure(1)
    for resolution in results:
        if not resolution[4:] in resolutions:
            resolutions.append(resolution[4:])

        log_z = []
        proc_time = []
        per_step_time = []
        for _, run in results[resolution].items():
            log_z.append(run["log_estimated_evidence"])
            proc_time.append(run["nestor_process_time"])
            per_step_time.append(run["mcmc_step_time"])

        all_log_z[resolution] = log_z
        mean_proc_time.append(np.mean(proc_time))
        mean_per_step_time.append(np.mean(per_step_time))

        avg_logz = np.mean(log_z)
        stderr_logz = np.std(log_z) / math.sqrt(len(log_z))
        plt.errorbar(
            resolution[4:], avg_logz, yerr=stderr_logz, fmt="o", c="dodgerblue"
        )

    plt.xlabel("Resolutions")
    plt.ylabel("log(Evidence)")
    plt.savefig(
        os.path.join(
            h_params["parent_path"],
            f"trial_{h_params['trial_name']}_evidence_errorbarplot.png",
        )
    )

    plt.figure(2)
    # resolutions, mean_proc_time = zip(
    #     *sorted(zip(resolutions, mean_proc_time), key=lambda x: x[0])
    # )
    plt.scatter(resolutions, mean_proc_time, c="C2", marker="o")
    plt.xlabel("Resolutions")
    plt.ylabel("Nested sampling process time")
    plt.savefig(
        os.path.join(
            h_params["parent_path"], f"trial_{h_params['trial_name']}_proctime.png"
        )
    )

    plt.figure(3)
    # resolutions, mean_per_step_time = zip(
    #     *sorted(zip(resolutions, mean_per_step_time), key=lambda x: x[0])
    # )
    plt.scatter(resolutions, mean_per_step_time, c="C2", marker="o")
    plt.xlabel("Resolutions")
    plt.ylabel("Mean time per MCMC step")
    plt.savefig(
        os.path.join(
            h_params["parent_path"], f"trial_{h_params['trial_name']}_persteptime.png"
        )
    )


def main(h_param_file, topology=True):
    with open(h_param_file, "r") as paramf:
        h_params = yaml.safe_load(paramf)

    max_allowed_runs = h_params["max_usable_threads"] // h_params["num_cores"]
    parent_path = h_params["parent_dir"]

    if not os.path.isdir(parent_path):
        os.mkdir(parent_path)

    target_runs = str(h_params["num_runs"])
    if "-" not in target_runs:
        target_runs = range(0, int(target_runs))
    else:
        target_runs = range(
            int(target_runs.split("-")[0]), int(target_runs.split("-")[1])
        )

    torun = get_all_toruns(h_params, target_runs)

    results = {}

    processes = {"Dummy": "Dummy"}
    completed_runs = []
    while len(list(processes.keys())) > 0:
        if "Dummy" in processes.keys():
            processes.pop("Dummy")

        if len(torun) > 0:
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
                        run_cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                    processes[(res, run_id)] = p
                    torun.remove((res, run_id))
                    print(f"Launched: {res.split('/')[-1]}, run_{run_id}")

                else:
                    print("Waiting for free threads...")

        waiting = True
        while waiting:
            for _, p in processes.items():
                if p.poll() is not None:
                    waiting = False
                # print(p.poll())

        (
            processes,
            curr_faulty_runs,
            successful_runs,
        ) = get_curr_processes_and_terminated_runs(processes)

        for proc in successful_runs:
            completed_runs.append(proc)
        if len(processes) == 0:
            break

        if len(curr_faulty_runs) != 0:
            for fr, p in curr_faulty_runs:
                if p.returncode == 11:
                    print(f"Will relaunch ({fr[0].split('/')[-1]}, run_{fr[1]})")
                    torun.append(fr)
                elif p.returncode == 12:
                    print(
                        f"Terminated: {fr[0].split('/')[-1]}, run_{fr[1]} with exit code: {p.returncode}"
                    )
                    print(
                        f"{fr[0].split('/')[-1]}, run_{fr[1]} ran out of maximum allowed iterations before converging. Will not relaunch it..."
                    )

    print(f"Waiting for {len(processes.keys())} processes to terminate...")

    while len(processes) > 0:
        final_waiting = True
        while final_waiting:
            for _, p in processes.items():
                if p.poll() is not None:
                    final_waiting = False

        (
            processes,
            curr_faulty_runs,
            successful_runs,
        ) = get_curr_processes_and_terminated_runs(processes)

        for proc in successful_runs:
            completed_runs.append(proc)

    ############## Preparing the output ###############

    print("Performing housekeeping tasks")

    for proc in completed_runs:
        run_deets, p = proc
        if p.returncode == 0:
            out, _ = p.communicate()

            result = literal_eval(out[4:])

            if run_deets[0].split("/")[-1] not in results.keys():
                results[f"{run_deets[0].split('/')[-1]}"] = {
                    f"run_{run_deets[1]}": result
                }
            else:
                results[f"{run_deets[0].split('/')[-1]}"][
                    f"run_{run_deets[1]}"
                ] = result

        else:
            _, err = p.communicate()
            print(err)
            exit()

    if "nestor_output.yaml" in os.listdir(parent_path):
        with open(os.path.join(parent_path, "nestor_output.yaml"), "r") as inf:
            old_results = yaml.safe_load(inf)
            merge(results, old_results)

    with open(f"{parent_path}/nestor_output.yaml", "w") as outf:
        yaml.dump(results, outf)


###################################################
###################### Main #######################
###################################################

if __name__ == "__main__":
    h_param_file = sys.argv[1]

    if sys.argv[2] == "manual":
        use_topology = False
    elif sys.argv[2] == "topology":
        use_topology = True

    with open(h_param_file, "r") as h_paramf:
        h_params = yaml.safe_load(h_paramf)
    if sys.argv[3] != "skip_calc":
        main(h_param_file, use_topology)

    else:

        with open(
            os.path.join(h_params["parent_path"], "nestor_output.yaml"), "r"
        ) as outf:
            results = yaml.safe_load(outf)

    if len(list(results.keys())) > 0:
        plotter(results, h_params)
    else:
        print("\nNone of the runs was successful...!")
    print("Done...!\n\n")
