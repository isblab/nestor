import os
import sys
import yaml
import shutil
import subprocess

h_param_file = sys.argv[1]
parent_dir = sys.argv[2]

with open(h_param_file,'r') as paramf:
    h_params = yaml.safe_load(paramf)

faulty_resolutions = []
faulty_runs = []

while 'faulty_runs.log' in os.listdir(parent_dir):
    with open(f'{parent_dir}/faulty_runs.log','r') as frf:
        for ln in frf.readlines():
            faulty_runs.append(ln.strip())
            faulty_resolutions.append(ln.strip().split('/')[-2])


    ##### Handle IMP runs
    processes = []
    for i,res in enumerate(faulty_resolutions):
        runid = faulty_runs[i].split('/')[-1]
        os.chdir(f'{parent_dir}/{res}')
        shutil.rmtree(runid)
        print(runid)
        os.system(f"mkdir run_{runid.split('_')[-1]}")
        os.chdir(f"run_{runid.split('_')[-1]}")

        run_cmd = ['mpirun','-n','4',
                                h_params['imp_path'], 'python', h_params['modeling_script_path'],
                                str(runid.split('_')[-1]), res.split('_')[-1], h_param_file]
        p = subprocess.Popen(run_cmd, stderr=subprocess.PIPE)
        processes.append(p)

    with open(f'{parent_dir}/faulty_runs_error.txt','w') as errf:
        for proc in processes:
            op,er = proc.communicate()
            er = er.decode()

            er2write = ''
            for erline in er.split('\n'):
                if not 'warning' in erline.lower():
                    er2write = er2write+'\n'+erline

            errf.write(str(proc.args))
            errf.write(er2write+'\n\n')

    os.remove(f'{parent_dir}/faulty_runs.log')
    subprocess.run(['python3','/home/shreyas/Projects/cgopt/optimising_rep/test_run_completion.py',h_param_file,parent_dir])
