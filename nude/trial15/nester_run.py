import os,sys
import subprocess

resolutions = ['01','05','10','20','30','50']
num_runs = 5

for res in resolutions:
    topology_file = f"topology{str(res)}.txt"
    os.system(f'mkdir res_{res}')
    os.chdir(f'res_{res}')
    for runid in range(num_runs):
        os.system(f'mkdir run_{runid}')
        os.chdir(f'run_{runid}')
        subprocess.Popen(['/home/shreyas/Projects/cgopt/imp-clean2/build/setup_environment.sh', 'python', '../../nude_modeling.py',str(runid),topology_file])
        os.chdir('../')
    os.chdir('../')
