import os
import sys
import glob
import math
import numpy as np
import matplotlib.pyplot as plt

def get_evidences_H(resolution_dir):
    run_dirs = glob.glob(os.path.join(resolution_dir,'run*'))
    evidences = []
    analytical_uncertainties = []
    for run in run_dirs:
        run_log_file = os.path.join(run,'run.log')
        with open(run_log_file,'r') as rlf:
            for ln in rlf.readlines():
                if ln.startswith('Accumulated evidence'):
                    evidences.append(float(ln.strip().split(': ')[-1]))
                if ln.startswith('Analytical uncertainty'):
                    analytical_uncertainties.append(float(ln.strip().split(': ')[-1]))
    return evidences, analytical_uncertainties


resolutions = sys.argv[1:]

x_vals = []
ana_unc_vals = []
evi_std_vals = []
for res in resolutions:
    evidences, ana_unc = get_evidences_H('res_'+res+'/')
    evidence_std = np.std(evidences)
    ana_unc_mean = np.mean(ana_unc)

    x_vals.append(f'res_{res}')
    ana_unc_vals.append(ana_unc_mean)
    evi_std_vals.append(math.log(evidence_std))

fig,ax1 = plt.subplots()
ax1.plot(x_vals,ana_unc_vals,c='C1',label='Analytical uncertainties')
ax1.set_xlabel('Resolutions')
ax1.set_ylabel('Analytical uncertainties')

ax2 = ax1.twinx()
ax2.plot(x_vals,evi_std_vals,c='C2',label='Log(std(evidences))')
ax2.set_ylabel('Log standard deviations of evidences')
fig.legend()
fig.savefig('tmp.png')
