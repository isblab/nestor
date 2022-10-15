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
            # print('Some error occurred in the run. Skipping that run...')

    return evidences, analytical_uncertainties


resolutions = sys.argv[1:]

x_vals = []
ana_unc_vals = []
ana_unc_mean_vals = []
evidence_mean_vals = []
evi_std_err_vals = []

for res in resolutions:
    evidences, ana_unc = get_evidences_H('res_'+res+'/')
    evidence_std_err = np.std(evidences)/math.sqrt(len(evidences))

    evidence_mean_vals.append(np.mean(evidences))

    for unc in ana_unc:
        x_vals.append(f'res_{res}')
        ana_unc_vals.append(unc)

    ana_unc_mean_vals.append(np.mean(ana_unc))
    evi_std_err_vals.append(evidence_std_err)
    # print(f"Analytical uncertainty standard deviation: {np.std(ana_unc)}") #/math.sqrt(len(ana_unc))}


uniq_x = []
for x in x_vals:
    if x not in uniq_x:
        uniq_x.append(x)


fig,ax = plt.subplots(3, sharex=True)

ax[0].errorbar(uniq_x, evidence_mean_vals, yerr=ana_unc_mean_vals, marker='o', label='logZ with ana_unc errorbar')
ax[0].set_ylabel('Log evidence')
# ax[0].legend()

ax[1].errorbar(uniq_x, evidence_mean_vals, yerr=evi_std_err_vals, marker='o', c='C3', label='logZ with std error errorbar')
ax[1].set_ylabel('Log evidence')
# ax[1].legend()

ax[2].scatter(x_vals,ana_unc_vals, c='C2', marker='o', label='Analytical uncertainties')
ax[2].plot(uniq_x, ana_unc_mean_vals, marker='*', c='C2', label='Mean analytical uncertainties')
ax[2].set_xlabel('Resolutions')
ax[2].set_ylabel('Analytical uncertainties')

ax2 = ax[2].twinx()
ax2.plot(uniq_x, evi_std_err_vals, c='C4', marker='o', label='Standard error (log(Z)')
ax2.set_ylabel('Stderr(log(evidences))')
# ax[2].legend()
# ax2.legend()

fig.legend()
fig.savefig('uncertainty_plot.png',dpi=400)
# plt.show()
