import os
import sys
import math
import yaml
import glob
import shutil
import linecache
import subprocess
import numpy as np
from matplotlib import pyplot as plt

###################################################
###################### Init #######################
###################################################

h_param_file = '/home/shreyas/Projects/cgopt/systems/to_ms/ninit/rnapolii/nestor_params_init_playground.yaml'
topology = False

if 'topology' in sys.argv:
	topology = True

with open(h_param_file,'r') as paramf:
	h_params = yaml.safe_load(paramf)


max_allowed_runs = h_params['max_usable_threads'] // h_params['num_cores']
imp_path = None
if 'imp_path' in h_params.keys():
    imp_path = h_params['imp_path']


###################################################
#################### Functions ####################
###################################################

def get_all_toruns(h_params):
	parent_dir = h_params['parent_dir']
	runs = []
	for res in h_params['resolutions']:
		for run in range(h_params['num_runs']):
			run_deets = (os.path.join(parent_dir,f"res_{res}"), str(run))
			runs.append(run_deets)
	return runs


def terminate_process(process):
	with open('error.txt','a') as errf:
		op, er = process.communicate()
		er = er.decode()
		er2write = ''
		for erline in er.split('\n'):
			if not 'warning' in erline.lower():
				if not erline.startswith('['):
					er2write = er2write+'\n'+erline
					errf.write(str(process.args))
					errf.write(er2write+'\n\n')


def fetch_process_deets_for_terminated_process(process):
	all_args = process.args
	res_id = os.path.join(h_params['parent_dir'], 'res_'+all_args[-2].split('.')[-2][-2:])
	run_id = process.args[-3]
	return res_id, run_id


def communicate_finished_proc_and_get_remaining_procs(processes, lastcall=False):
	terminated_processes = []
	faulty_runs = []
	if not lastcall:
		for proc in processes:
			if proc.poll() is not None:
				terminate_process(proc)
				processes.remove(proc)
				terminated_processes.append(proc.poll())
					
				if proc.poll()==11:
					termination_deets = fetch_process_deets_for_terminated_process(proc)
					faulty_runs.append((termination_deets, proc.poll()))
					shutil.rmtree(os.path.join(termination_deets[0],'run_'+termination_deets[1]))
					
	else:
		for proc in processes:
			terminate_process(proc)
			terminated_processes.append(proc.poll())
			if not proc.poll() is None:
				if proc.poll()!=0:
					faulty_runs.append((fetch_process_deets_for_terminated_process(proc), proc.poll()))
			processes.remove(proc)
	        
	return processes, terminated_processes, faulty_runs


def get_torun(all_runs):
	to_run = []
	for run_deets in all_runs:
		path_to_look = os.path.join(run_deets[0],f"run_{run_deets[1]}/")
		if not os.path.isdir(path_to_look):
			to_run.append(run_deets)
	return to_run


def concatenate_evidences(resolutions):
    c_ev_files = []
    for res in resolutions:
        dir_name = 'res_'+res
        os.chdir(dir_name)
        evidence_files = glob.glob('run*/estimated_evidence*')
        with open("estimated_evidences.dat",'w') as concatf:
            for fl in evidence_files:
                with open(fl,'r') as evf:
                    concatf.write(evf.read())
        os.chdir('../')
        c_ev_files.append(f"{dir_name}/estimated_evidences.dat")
    return c_ev_files


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


def plot_evidence_w_stderr(h_params):
    files = concatenate_evidences(h_params['resolutions'])
    all_log_evi = []
    output = {}

    for res in files:
        evidences = []
        log_evidences = []

        with open(res,'r') as evf:
            for ln in evf.readlines():
                evidence = float(ln.strip())
                evidences.append(evidence)
                log_evidences.append(math.log(evidence))
                all_log_evi.append(math.log(evidence))
        std_err = np.std(log_evidences)/math.sqrt(len(log_evidences))
        plt.errorbar(int(res.split('/')[0][-2:]), np.mean(log_evidences),
                    yerr=std_err, fmt='o')

        process_times = get_process_time(res.split('/')[-2])
        if 'nestOR_output.log' in os.listdir(h_params['parent_dir']):
            with open('nestOR_output.log','r') as outl:
                output = yaml.safe_load(outl)

        output[f"Resolution {res.split('/')[0][-2:]}"] = {'Mean log evidence':float(np.mean(log_evidences)),
                                                            'Standard error on log evidence':float(std_err),
                                                            'Time taken':float(np.mean(process_times))}
        with open('nestOR_output.log','w') as outl:
            output = yaml.dump(output, outl, default_flow_style = False)

    #plt.yticks(np.arange(int(min(all_log_evi)-5),int(max(all_log_evi))+5,2))
    # plt.grid(axis='y')
    plt.savefig(f"trial_{h_params['trial_name']}_evidence.png")


###################################################
###################### Main #######################
###################################################

terminate = False
all_runs = get_all_toruns(h_params)
torun = [run for run in all_runs]
print(torun)

max_allowed_runs = h_params['max_usable_threads'] // h_params['num_cores']
imp_path = None
if 'imp_path' in h_params.keys():
    imp_path = h_params['imp_path']


processes = []
while len(torun) != 0:
    resolutions = [ids[0] for ids in torun]
    for i, res in enumerate(resolutions):
        # Run all processes that can be run on the current number of cores
        os.chdir(h_params['parent_dir'])
        if not os.path.isdir(res):
            os.system(f'mkdir {res}')

        if len(processes) < max_allowed_runs:
            if topology:
                topf = f"topology{res.split('/')[-1].split('_')[-1]}.txt"
            else:
                topf = res.split('/')[-1].split('_')[-1]

            os.chdir(res)
            if imp_path is None:
                run_cmd = ['mpirun','-n',str(h_params['num_cores']),
                            'python', h_params['modeling_script_path'],
                            str(torun[i][1]), topf, h_param_file]
            else:
                run_cmd = ['mpirun','-n',str(h_params['num_cores']), h_params['imp_path'],
                            'python', h_params['modeling_script_path'],
                            str(torun[i][1]), topf, h_param_file]

            os.system(f'mkdir run_{torun[i][1]}')
            os.chdir(f'run_{torun[i][1]}')
            
            p = subprocess.Popen(run_cmd, stderr=subprocess.PIPE)
            processes.append(p)
            os.chdir('../../')
        
    terminated_processes = []
    while len(terminated_processes)==0:
    	processes, terminated_processes, _ = communicate_finished_proc_and_get_remaining_procs(processes)
    	
    torun = get_torun(all_runs)
    with open('tr.txt','a') as tr:
    	for t in torun:
    		tr.write('An iteration completed ' + t[0] + ',' + t[1])
    	tr.write('\n\n')
print('Finishing up')
processes, terminated_processes, faulty_runs = communicate_finished_proc_and_get_remaining_procs(processes, lastcall=True)
plot_evidence_w_stderr(h_params)
print('Done with the runs')