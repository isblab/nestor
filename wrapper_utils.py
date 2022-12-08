import os,sys
import yaml
import shutil
import linecache


def generate_initial_torun():
    torun = []
    nruns = h_params['num_runs']
    parent_dir = h_params['parent_dir']
    resolutions = h_params['resolutions']
    for res in resolutions:
        for run in range(nruns):
            torun.append((os.path.join(parent_dir,f"res_{res}"), run))
    return torun


def test_run_completion(all_runs):
    to_be_run = []
    runs = []

    for i in all_runs:
        runs.append(os.path.join(str(i[0]),'run_'+str(i[1])))

    for run in runs:
        if not os.path.isdir(run):
            out_run = ('/'.join(run.split('/')[:-1]), int(run.split('/')[-1].split('_')[-1]))
            if not out_run in to_be_run:
                to_be_run.append(out_run)

        if os.path.isdir(run) and 'error.log' in os.listdir(run):
            error = linecache.getline(f'{run}/error.log',3).split(':')[-1].strip()
            if not 'MaxIterations reached without convergence criteria' in error:
                out_run = ('/'.join(run.split('/')[:-1]), int(run.split('/')[-1].split('_')[-1]))
                # if not out_run in to_be_run:
                to_be_run.append(out_run)
                shutil.rmtree(f"{out_run[0]}/run_{out_run[1]}")

    return to_be_run
