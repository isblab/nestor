import os
import shutil
import IMP.nestor.wrapper_v6 as wrapper_v6

expected_files = [
    "trial_optrep_params_evidence_errorbarplot.png",
    "trial_optrep_params_persteptime.png",
    "trial_optrep_params_proctime.png",
    "nestor_output.yaml",
]

if "runs" in os.listdir(os.getcwd()):
    shutil.rmtree(os.path.join(os.getcwd(), "runs"))

paramf_path = "test/input/nestor_params_optrep.yaml"
wrapper_exitcode = wrapper_v6.main(paramf_path, True)
