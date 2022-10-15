import os
import sys
import glob
import math
import yaml
import numpy as np
import matplotlib.pyplot as plt


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
            pass

    return evidences, analytical_uncertainties


runs = sys.argv[1:]
resolutions = ['01','05','10','20','30','50']
fig, ax = plt.subplots(2, sharex=True)

master_stderr_evidences = []
master_stderr_ana_uncertainties = []

for run in runs:
    xvals = []
    stderr_evidences = []
    stderr_ana_uncertainties = []

    for res in resolutions:
        evidences, ana_unc = get_evidences_H(run+'/res_'+res+'/')
        stderr_evi = np.std(evidences) / (math.sqrt(len(evidences)))
        stderr_evidences.append(stderr_evi)
        stderr_ana_unc = np.std(ana_unc) / (math.sqrt(len(ana_unc)))
        stderr_ana_uncertainties.append(stderr_ana_unc)
        xvals.append(f'res_{res}')

    master_stderr_evidences.append(stderr_evidences)
    master_stderr_ana_uncertainties.append(stderr_ana_uncertainties)



sterrevi, = ax[0].plot(xvals, master_stderr_evidences[0], marker='o', label=runs[0], c='C1')
ax[0].plot(xvals, master_stderr_evidences[1], marker='o', label=runs[1], c='C2')

ax[1].plot(xvals, master_stderr_ana_uncertainties[0], marker='o', label=runs[0], c='C1')
sterrana, = ax[1].plot(xvals, master_stderr_ana_uncertainties[1], marker='o', label=runs[1], c='C2')

ax[0].set_ylabel('StdErr(log(Z))')
ax[1].set_xlabel('Resolutions')
ax[1].set_ylabel('StdErr(analytical uncertainty)')


fig.legend(handles=[sterrevi,sterrana])
fig.savefig(f'comparision between {runs[0]} and {runs[1]}.png',dpi=600)
