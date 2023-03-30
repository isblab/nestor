import os
import yaml
import pandas as pd

parent_dir = "/home/shreyas/Projects/cgopt/optimising_rep/utils/benchmarks/"
fnames = (
    os.path.join(parent_dir, "gtusc.csv"),
    os.path.join(parent_dir, "rnapolii.csv"),
    os.path.join(parent_dir, "mhm.csv"),
    os.path.join(parent_dir, "nude.csv"),
)

out_results = {}
for fname in fnames:
    results = {}
    df = pd.read_csv(fname)
    for i in range(df.shape[0]):
        benchmark_deets = df.loc[i].to_dict()
        resolution = benchmark_deets["Resolution"]
        del benchmark_deets["Resolution"]
        results[resolution] = benchmark_deets
    out_results[fname.split("/")[-1].split(".")[-2]] = results

with open(os.path.join(parent_dir, "benchmarks_combined.yaml"), "w") as outf:
    yaml.dump(out_results, outf)
