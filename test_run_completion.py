import os
import sys
import yaml
import glob
import linecache


h_param_file = sys.argv[1]      # param file
parent_dir = sys.argv[2]        # directory that contains all resolutions

with open(h_param_file,'r') as paramf:
    h_params = yaml.safe_load(paramf)

faulty_runs = []
resolution_dirs_created = glob.glob(f'{parent_dir}/res*')

for res in resolution_dirs_created:
    runs = glob.glob(f'{res}/run*')
    for run in runs:
        if 'error.log' in os.listdir(run):
            error = linecache.getline(f'{run}/error.log',3).split(':')[-1].strip()
            if not 'MaxIterations reached without convergence criteria' in error:
                faulty_runs.append(run)

if len(faulty_runs)!=0:
    with open(f'{parent_dir}/faulty_runs.log','w') as frf:
        for runid in faulty_runs:
            frf.write(f'{runid}\n')
