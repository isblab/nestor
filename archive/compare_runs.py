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


runs = sys.argv[1:]
resolutions = ['01','05','10','20','30','50']
fig, ax = plt.subplots(3, sharex=True)

fig1 = plt.figure()
fax1 = fig1.add_subplot()
fig2 = plt.figure()
fax2 = fig2.add_subplot()
fig3 = plt.figure()
fax3 = fig3.add_subplot()
fig4 = plt.figure()
fax4 = fig4.add_subplot()


master_stderr_evidences = []
master_stderr_ana_uncertainties = []
t_plots = []

for i,run in enumerate(runs):
    xvals = []
    times = []
    mean_evidences = []
    stderr_evidences = []
    stderr_ana_uncertainties = []

    for res in resolutions:
        times.append(np.mean(get_process_time(run+'/res_'+res+'/')))
        evidences, ana_unc = get_evidences_H(run+'/res_'+res+'/')
        stderr_evi = np.std(evidences) / (math.sqrt(len(evidences)))
        stderr_evidences.append(stderr_evi)
        stderr_ana_unc = np.std(ana_unc) / (math.sqrt(len(ana_unc)))
        stderr_ana_uncertainties.append(stderr_ana_unc)
        xvals.append(f'res_{res}')
        mean_evidences.append(np.mean(evidences))

    master_stderr_evidences.append(stderr_evidences)
    master_stderr_ana_uncertainties.append(stderr_ana_uncertainties)

    fax3.plot(xvals,times, marker='o',label=run, c=f'C{i+1}')
    fax4.errorbar(xvals,mean_evidences,yerr=stderr_evidences, marker='o',label=run, c=f'C{i+1}')


for i,run in enumerate(runs):
    fax1.plot(xvals, master_stderr_evidences[i], marker='o', label=runs[i], c=f'C{i+1}')
# fax1.plot(xvals, master_stderr_evidences[1], marker='o', label=runs[1], c='C2')
# fax1.plot(xvals, master_stderr_evidences[2], marker='o', label=runs[2], c='C3')

    fax2.plot(xvals, master_stderr_ana_uncertainties[i], marker='o', label=runs[i], c=f'C{i+1}')
# fax2.plot(xvals, master_stderr_ana_uncertainties[1], marker='o', label=runs[1], c='C2')
# fax2.plot(xvals, master_stderr_ana_uncertainties[2], marker='o', label=runs[2], c='C3')

fax1.set_xlabel('Resolutions')
fax1.set_ylabel('StdErr(log(Z))')

fax2.set_xlabel('Resolutions')
fax2.set_ylabel('StdErr(analytical uncertainty)')

fax3.set_xlabel('Resolutions')
fax3.set_ylabel('Time taken (seconds)')

fax4.set_xlabel('Resolutions')
fax4.set_ylabel('Evidences with std_err errorbar')

fax1.legend()
fax2.legend()
fax3.legend()
fax4.legend()

fig1.savefig(f'stderr_comparision.png',dpi=600)
fig2.savefig(f'anaerr_comparision.png',dpi=600)
fig3.savefig(f'timetaken_comparision.png',dpi=600)
fig4.savefig(f'evidences_comparision.png',dpi=600)
