import os,sys
import time
import glob
import math
import yaml
import shutil
import linecache
import subprocess
import numpy as np
from matplotlib import pyplot as plt

topology = False

h_param_file = sys.argv[1]
if 'topology' in sys.argv:
    topology = True


# TODO borrow from other scripts, make it as functions in this script.
# TODO while runs_reqd is not empty, do runs. test runs
# TODO runs_reqd = [(runid,res),...]
# TODO minimize file write /output


# -----------------------------------------------------------------------------
# --------------------------------- Functions ---------------------------------
# -----------------------------------------------------------------------------

def generate_initial_torun():
    torun = []
    nruns = h_params['num_runs']
    parent_dir = h_params['parent_dir']
    resolutions = h_params['resolutions']
    for res in resolutions:
        for run in range(nruns):
            torun.append((os.path.join(parent_dir,f"res_{res}"), run))
    return torun


def get_finished_pid(processes):
    finished_process_id = []
    for idx, proc in enumerate(processes):
        if proc.poll() is not None:
            finished_process_id.append(idx)
    return finished_process_id


def get_remaining_processes(processes,finish=False):
    remaining_processes = processes.copy()
    with open('error.txt','a') as errf:
        if finish:
            finished_process_ids = range(len(processes))
        else:    
            finished_process_ids = get_finished_pid(processes)

        for pid in finished_process_ids:
            proc = processes[pid]
            remaining_processes.remove(proc)

            op,er = proc.communicate()
            er = er.decode()
            er2write = ''
            for erline in er.split('\n'):
                if not 'warning' in erline.lower():
                    if not erline.startswith('['):
                        er2write = er2write+'\n'+erline
            errf.write(str(proc.args))
            errf.write(er2write+'\n\n')

    return remaining_processes


def test_run_completion(resolutions):
    all_runs = generate_initial_torun()
    faulty_runs = []

    runs = []
    for i in all_runs:
        runs.append(os.path.join(str(i[0]),'run_'+str(i[1])))

    for run in runs:
        if not os.path.isdir(run):
            out_run = ('/'.join(run.split('/')[:-1]), int(run.split('/')[-1].split('_')[-1]))
            if not out_run in faulty_runs:
                faulty_runs.append(out_run)

        elif 'error.log' in os.listdir(run):
            error = linecache.getline(f'{run}/error.log',3).split(':')[-1].strip()
            if not 'MaxIterations reached without convergence criteria' in error:
                out_run = ('/'.join(run.split('/')[:-1]), int(run.split('/')[-1].split('_')[-1]))
                if not out_run in faulty_runs:
                    faulty_runs.append(out_run)
                    shutil.rmtree(f"{out_run[0]}/run_{out_run[1]}")


    return faulty_runs


def concatenate_evidences(resolutions):
    c_ev_files = []
    for res in resolutions:
        dir_name = 'res_'+res
        os.chdir(dir_name)
        evidence_files = glob.glob('run*/estimated_evidence*')
        with open("estimated_evidences.dat",'w') as concatf:
            for fl in evidence_files:
                with open(fl,'r') as evf:
                    concatf.write(evf.read())
        os.chdir('../')
        c_ev_files.append(f"{dir_name}/estimated_evidences.dat")
    return c_ev_files


def get_process_time(resolution_dir):
    run_dirs = glob.glob(os.path.join(resolution_dir,'run*'))
    times = []
    for run in run_dirs:
        run_log_file = os.path.join(run,'run.log')
        try:
            with open(run_log_file,'r') as rlf:
                for ln in rlf.readlines():
                    if ln.startswith('Nestor process time:'):
                        times.append(float(ln.strip().split(': ')[-1].split(' ')[0]))
        except FileNotFoundError:
            print('Found a directory with no log file')

    return times


def plot_evidence_w_stderr(h_params):
    files = concatenate_evidences(h_params['resolutions'])
    all_log_evi = []
    output = {}

    for res in files:
        evidences = []
        log_evidences = []

        with open(res,'r') as evf:
            for ln in evf.readlines():
                evidence = float(ln.strip())
                evidences.append(evidence)
                log_evidences.append(-math.log(evidence))
                all_log_evi.append(-math.log(evidence))
        std_err = np.std(log_evidences)/math.sqrt(len(log_evidences))
        plt.errorbar(int(res.split('/')[0][-2:]), np.mean(log_evidences),
                    yerr=std_err, fmt='o')

        process_times = get_process_time(res.split('/')[-2])
        if 'nestOR_output.log' in os.listdir(h_params['parent_dir']):
            with open('nestOR_output.log','r') as outl:
                output = yaml.safe_load(outl)

        output[f"Resolution {res.split('/')[0][-2:]}"] = {'Mean log evidence':float(np.mean(log_evidences)),
                                                            'Standard error on log evidence':float(std_err),
                                                            'Time taken':float(np.mean(process_times))}
        with open('nestOR_output.log','w') as outl:
            output = yaml.dump(output, outl, default_flow_style = False)

    #plt.yticks(np.arange(int(min(all_log_evi)-5),int(max(all_log_evi))+5,2))
    # plt.grid(axis='y')
    plt.savefig(f"trial_{h_params['trial_name']}_evidence.png")


# -----------------------------------------------------------------------------
# ----------------------------------- Main ------------------------------------
# -----------------------------------------------------------------------------

# ----- Init -----

with open(h_param_file, 'r') as paramf:
    h_params = yaml.safe_load(paramf)

torun = generate_initial_torun()

max_allowed_runs = h_params['max_usable_threads'] // h_params['num_cores']
imp_path = None
if 'imp_path' in h_params.keys():
    imp_path = h_params['imp_path']

# ----- IMP runs -----

processes = []

while len(torun) != 0:
    resolutions = [ids[0] for ids in torun]
    for i, res in enumerate(resolutions):
        os.chdir(h_params['parent_dir'])
        if not res.split('/')[-1] in os.listdir(h_params['parent_dir']):
            os.system(f'mkdir {res}')

        if len(processes) < max_allowed_runs:
            with open('ran.txt','a') as rt:
                rt.write(f"{torun[i]}\n")
            if topology:
                topf = f"topology{res.split('/')[-1].split('_')[-1]}.txt"
            else:
                topf = res.split('/')[-1].split('_')[-1]

            os.chdir(res)
            if imp_path is None:
                run_cmd = ['mpirun','-n',str(h_params['num_cores']),
                            'python', h_params['modeling_script_path'],
                            str(torun[i][1]), topf, h_param_file]
            else:
                run_cmd = ['mpirun','-n',str(h_params['num_cores']), h_params['imp_path'],
                            'python', h_params['modeling_script_path'],
                            str(torun[i][1]), topf, h_param_file]
                
            with open(f"{h_params['parent_dir']}/runs.txt",'a') as rtf:
                rtf.write(f"{torun[i]}\n")
            if 'run_{torun[i][1]}' in os.listdir('./'):
                shutil.rmtree('run_{torun[i][1]}')
                with open(f"{h_params['parent_dir']}/runs.txt",'a') as rtf:
                    rtf.write(f"\nSpawning new process: {run_cmd}\n\n")

            os.system(f'mkdir run_{torun[i][1]}')
            os.chdir(f'run_{torun[i][1]}')
            p = subprocess.Popen(run_cmd, stderr=subprocess.PIPE)
            processes.append(p)
            os.chdir('../../')
        else:
            with open('ran.txt','a') as rt:
                rt.write(f"\n")

    while len(get_finished_pid(processes))==0:
        pass
    processes = get_remaining_processes(processes)
    torun = test_run_completion(resolutions)
    
get_remaining_processes(processes,finish=True)
plot_evidence_w_stderr(h_params)
print('Done with the runs')
