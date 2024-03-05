import os
import shutil
import subprocess

expected_files = [
    "trial_optrep_params_evidence_errorbarplot.png",
    "trial_optrep_params_persteptime.png",
    "trial_optrep_params_proctime.png",
    "nestor_output.yaml",
]

if "runs" in os.listdir(os.getcwd()):
    shutil.rmtree(os.path.join(os.getcwd(), "runs"))

command = [
    "python",
    "../pyext/src/wrapper_v6.py",
    "input/nestor_params_optrep.yaml",
]
p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
if p.returncode != 0:
    print(err)
    print(out)

generated_files = os.listdir(os.path.join(os.getcwd(), "runs"))

for exp_file in expected_files:
    print(exp_file in generated_files)
