import os,sys
import math
import numpy as np
from matplotlib import pyplot as plt

trial_id = sys.argv[1]
files = ['res_01/estimated_evidences.dat',\
        'res_05/estimated_evidences.dat',\
        'res_10/estimated_evidences.dat',\
        'res_20/estimated_evidences.dat',\
        'res_30/estimated_evidences.dat',\
        'res_50/estimated_evidences.dat']

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
    plt.errorbar(res.split('/')[0], np.mean(log_evidences), yerr=std_err, fmt='o')

plt.yticks(np.arange(int(min(mean_log_evi)-10),int(max(mean_log_evi))+10,2))
plt.grid(axis='y')
plt.savefig(f'trial_{trial_id}_evidence.png') # Include trial number in fname
plt.show()
