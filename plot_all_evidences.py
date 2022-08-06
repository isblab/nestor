import os
import sys
import glob
import numpy as np
from math import log, sqrt
from matplotlib import pyplot as plt

def get_evidences(directory):
    all_res_dir = glob.glob(os.path.join(directory,'res*'))
    all_evidences = {}

    for res in all_res_dir:
        with open(os.path.join(res,'estimated_evidences.dat'),'r') as evf:
            evidences = []
            for ln in evf.readlines():
                evidences.append(log(float(ln.strip())))

            all_evidences[res.split('/')[-1]] = evidences
    return all_evidences


trials = sys.argv[1].split(',')
for trial in trials:
    evidences = get_evidences(trial)
    std_err = []
    means = []
    for k in evidences:
        std_err.append(np.std(evidences[k])/sqrt(len(evidences[k])))
        means.append(np.mean(evidences[k]))

    plt.errorbar(evidences.keys(), means, yerr=std_err, fmt='o-',label=trial.split('/')[-1])

plt.legend()
plt.show()
