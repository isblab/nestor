import os
import glob

conditions = ["ninit", "nfpi", "prior_robustness"]
system_rename = {"gtusc": "γTuSC", "rnapolii": "RNA polymerase II"}
parent_path = "/home/shreyas/Projects/cgopt/systems/to_ms/latest"


def remove_input_dirs(dir_lst: list[str]):
    out_dirs = []
    for i in dir_lst:
        if not "data" in i:
            if not "input" in i:
                out_dirs.append(i)
    return out_dirs


for condn in conditions:
    p_path = os.path.join(parent_path, condn, "*")
    systems = glob.glob(p_path)

    print(f"Working on {condn}")
    for system in systems:
        sys_name = system.split("/")[-1]
        print(f"\tWorking on {sys_name}")
        if sys_name in system_rename:
            ori_sysname = sys_name
            sys_name = system_rename[sys_name]
            print(f"\t\tRenaming {ori_sysname} as {sys_name}")
        runsets = [i for i in glob.glob(os.path.join(system, "*")) if os.path.isdir(i)]
        runsets = sorted(remove_input_dirs(runsets))
        cmd = f"python ~/Projects/cgopt/optimising_rep/utils/compare_runs_v2_w_pyplot.py '{sys_name}' "
        for runset in runsets:
            cmd += f"{runset.split('/')[-1]} "

        os.chdir(system)
        os.system(cmd)
