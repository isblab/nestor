import os
import sys
import glob
import math
import subprocess
import numpy as np
from matplotlib import pyplot as plt

# -----------------------------------------------------------------------------
# -------------------------------- User Inputs --------------------------------
# -----------------------------------------------------------------------------

# '01','05','10',
resolutions = [20,30,50]
num_runs = 3
trial_id = 'new_world'

imp_path = '/home/shreyasarvindekar/Projects/cgopt/imp-clean2/build/setup_environment.sh'
modeling_script_path = '../../nude_modeling.py'


# -----------------------------------------------------------------------------
# --------------------------------- Functions ---------------------------------
# -----------------------------------------------------------------------------
def concatenate_evidences(resolutions):
    c_ev_files = []
    for res in resolutions:
        dir_name = 'res_'+str(res)
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
processes = []
for res in resolutions:
    topology_file = f"topology{str(res)}.txt"
    os.system(f'mkdir res_{res}')
    os.chdir(f'res_{res}')
    for runid in range(num_runs):
        os.system(f'mkdir run_{runid}')
        os.chdir(f'run_{runid}')
        p = subprocess.Popen([imp_path,'python',modeling_script_path,
                                str(runid),topology_file],
                                stderr=subprocess.PIPE)
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

print('Done with the runs')


# --------------------------------- Plotting ----------------------------------

files = concatenate_evidences(resolutions)
mean_log_evi = []
for res in files:
    evidences = []
    log_evidences = []

    with open(res,'r') as evf:
        for ln in evf.readlines():
            evidence = float(ln.strip())
            evidences.append(evidence)
            log_evidences.append(-math.log(evidence))

    std_err = np.std(log_evidences)/math.sqrt(len(log_evidences))
    mean_log_evi.append(np.mean(log_evidences))
    plt.errorbar(int(res.split('/')[0][-2:]), np.mean(log_evidences),
                yerr=std_err, fmt='o')

plt.yticks(np.arange(int(min(mean_log_evi)-10),int(max(mean_log_evi))+10,2))
plt.grid(axis='y')
plt.savefig(f'trial_{trial_id}_evidence.png') # Include trial number in fname
plt.show()
