import os
import sys
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


def write_error_file(processes):
    with open('error.txt','a') as errf:
        for proc in processes:
            op,er = proc.communicate()
            er = er.decode()
            er2write = ''
            for erline in er.split('\n'):
                if not 'warning' in erline.lower():
                    if not erline.startswith('['):
                        er2write = er2write+'\n'+erline
            errf.write(str(proc.args))
            errf.write(er2write+'\n\n')
    return []


def test_run_completion(resolutions):
    faulty_runs = []
    for res in resolutions:
        runs = glob.glob(f'{res}/run*')
        for run in runs:
            if 'error.log' in os.listdir(run):
                error = linecache.getline(f'{run}/error.log',3).split(':')[-1].strip()
                if not 'MaxIterations reached without convergence criteria' in error:
                    out_run = ('/'.join(run.split('/')[:-1]), run.split('/')[-1])
                    faulty_runs.append(run)
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


# -----------------------------------------------------------------------------
# ----------------------------------- Main ------------------------------------
# -----------------------------------------------------------------------------

# ----- Init -----

with open(h_param_file, 'r') as paramf:
    h_params = yaml.safe_load(paramf)

torun = generate_initial_torun()


# ----- IMP runs -----

processes = []

while len(torun) != 0:
    resolutions = [ids[0] for ids in torun]
    for i, res in enumerate(resolutions):
        if not res.split('/')[-1] in os.listdir(h_params['parent_dir']):
            os.system(f'mkdir {res}')

        if len(processes) < h_params['max_allowed_processes'] - h_params['num_cores']:
            if topology:
                topf = f'topology{res}.txt'
            else:
                topf = res

            os.chdir(res)
            run_cmd = ['mpirun','-n',str(h_params['num_cores']), h_params['imp_path'],
                        'python', h_params['modeling_script_path'],
                        str(torun[i][1]), topf, h_param_file]

            if 'run_{torun[i][1]}' in os.listdir('./'):
                shutil.rmtree('run_{torun[i][1]}')

            os.system(f'mkdir run_{torun[i][1]}')
            os.chdir(f'run_{torun[i][1]}')
            p = subprocess.Popen(run_cmd, stderr=subprocess.PIPE)
            processes.append(p)
            os.chdir('../../')
        else:
            processes = write_error_file(processes)

    processes = write_error_file(processes)
    torun = test_run_completion(resolutions)

print('Done with the runs')


# --------------------------------- Plotting ----------------------------------

files = concatenate_evidences(h_params['resolutions'])
all_log_evi = []
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

#plt.yticks(np.arange(int(min(all_log_evi)-5),int(max(all_log_evi))+5,2))
# plt.grid(axis='y')
plt.savefig(f"trial_{h_params['trial_name']}_evidence.png")
