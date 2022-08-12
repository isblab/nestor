'''
Given an CSV file for XLs in Protein1,Residue1,Protein2,Residue2 format,
this script will split it in the given ratio for use in Nested Sampling
'''
import os,sys
import random

xl_file = sys.argv[1]
perc_to_evi = 0.7

xls = []
with open(xl_file,'r') as xlf:
    for ln in xlf.readlines():
        if (not ln.startswith('Protein1')) and (not ln.startswith('Linker')):
            xls.append(ln)
        else:
            header = ln

sampling, evi_calc = [], []
for link in xls:
    rng = random.random()
    if rng<perc_to_evi:
        evi_calc.append(link)
    else:
        sampling.append(link)

fname = xl_file.split('/')[-1]
dir_path = xl_file.split('/')
if len(dir_path)>1:
    dir_path = '/'.join(dir_path[0:-1])
else:
    dir_path = './'
with open(f'{dir_path}/sampling_{fname}','w') as sf:
    sf.write(header)
    for lnk in sampling:
        sf.write(lnk)

with open(f'{dir_path}/evicalc_{fname}','w') as evif:
    evif.write(header)
    for lnk in evi_calc:
        evif.write(lnk)
