import os,sys
import math
import glob
import numpy as np
from matplotlib import pyplot as plt

def get_evidences_H(resolution_dir):
    run_dirs = glob.glob(os.path.join(resolution_dir,'run*'))
    evidences = []
    analytical_uncertainties = []
    for run in run_dirs:
        run_log_file = os.path.join(run,'run.log')
        try:
            with open(run_log_file,'r') as rlf:
                for ln in rlf.readlines():
                    if ln.startswith('Accumulated evidence'):
                        evidences.append(math.log(float(ln.strip().split(': ')[-1])))
                    if ln.startswith('Analytical uncertainty'):
                        analytical_uncertainties.append(float(ln.strip().split(': ')[-1]))
        except FileNotFoundError:
            print('Shuffle configuration error found. Skipping that run...')

    return evidences, analytical_uncertainties


resolutions = sys.argv[1:]

x_vals = []
evi_mean_vals = []
evi_std_err_vals = []
for res in resolutions:
    try:
        evidences, _ = get_evidences_H('res_'+res+'/')
        evi_std_err = np.std(evidences)/math.sqrt(len(evidences))
        evi_mean = np.mean(evidences)

        x_vals.append(f'res_{res}')
        evi_mean_vals.append(evi_mean)
        evi_std_err_vals.append(evi_std_err)
        plt.errorbar(f'res_{res}', evi_mean, yerr=evi_std_err, marker='o')


    except FileNotFoundError:
        print('Shuffle configuration error found. Skipping that run...')

plt.show()
