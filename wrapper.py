import os
import sys
import glob
import math
import yaml
import subprocess
import numpy as np
from matplotlib import pyplot as plt


topology = False

h_param_file = sys.argv[1]
if 'topology' in sys.argv:
    topology = True


# -----------------------------------------------------------------------------
# --------------------------------- Functions ---------------------------------
# -----------------------------------------------------------------------------
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

with open(h_param_file, 'r') as paramf:
    h_params = yaml.safe_load(paramf)

if not 'resolution_sets' in h_params.keys():
    h_params['resolution_sets'] = (h_params['resolutions'],)


for res1 in h_params['resolution_sets']:
    processes = []
    for res in res1:
        if topology:
            topf = f'topology{res}.txt'
        else:
            topf = res
        os.system(f'mkdir res_{res}')
        os.chdir(f'res_{res}')
        for runid in range(h_params['num_runs']):
            run_cmd = ['mpirun','-n','4',
                            h_params['imp_path'], 'python', h_params['modeling_script_path'],
                            str(runid), topf, h_param_file]

            os.system(f'mkdir run_{runid}')
            os.chdir(f'run_{runid}')
            p = subprocess.Popen(run_cmd, stderr=subprocess.PIPE)
            processes.append(p)
            os.chdir('../')
        os.chdir('../')

    with open('error.txt','w') as errf:
        for proc in processes:
            op,er = proc.communicate()
            er = er.decode()

            er2write = ''
            for erline in er.split('\n'):
                if not 'warning' in erline.lower():
                    er2write = er2write+'\n'+erline

            errf.write(str(proc.args))
            errf.write(er2write+'\n\n')


# ---------------------------- Handle faulty runs -----------------------------

subprocess.run(['python3',
                '/home/shreyas/Projects/cgopt/optimising_rep/test_run_completion.py',
                h_param_file, './'])

subprocess.run(['python3',
                '/home/shreyas/Projects/cgopt/optimising_rep/rerun_faulty_runs.py',
                h_param_file, './'])

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
