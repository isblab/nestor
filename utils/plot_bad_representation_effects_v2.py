import os
import sys
import yaml
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

path2config = sys.argv[1]
with open(path2config, "r") as p2c:
    loaded = yaml.safe_load(p2c)
    rep_key = loaded["systems"]
    path2data = loaded["path2benchmark_data"]


def plot_bars(representation_key, datapath, field):
    for system, deets in representation_key.items():
        df = pd.read_csv(os.path.join(datapath, f"{system}.csv"))
        optrep_df = df.loc[df["Resolution"] == deets["Optimal representation"]]
        incorrect_rep_df = df.loc[df["Resolution"] == deets["Incorrect representation"]]

        try:
            opt_y = float(optrep_df[field])
            incorr_y = float(incorrect_rep_df[field])
        except KeyError:
            continue

        plt.figure(1, figsize=(4, 6))
        x1 = 0
        width = 0.2
        x2 = x1 + width
        plt.bar(x=x1, height=incorr_y, width=width, color="orange")
        plt.bar(x=x2, height=opt_y, width=width, color="green")
        plt.bar(x=x2 + 2 * width, height=0)

        plt.xticks([])
        plt.ylabel(field)
        plt.title(system)
        plt.xlim(-1 * width, x2 + 7 * width)
        # plt.box(False)
        # plt.axes().get_xaxis().set_visible(True)
        plt.savefig(
            f"Comparision of {field} between optrep and incorrect rep_{system}.png",
            dpi=1200,
        )
        plt.close()
        break


for field in ("Model Precision (Å)", "GMM CC", "Normalized XL score", "Time (sec)"):
    plot_bars(rep_key, path2data, field)
    break
