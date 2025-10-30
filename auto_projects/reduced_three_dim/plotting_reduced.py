import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_one_param(filename, figurename):
    data_eq = []
    data_periodic = []
    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            if parts and parts[0] == '1':
                try:
                    data_eq.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_eq.append(row)
            if parts and parts[0] == '2':
                try:
                    data_eq.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_eq.append(row)
            if parts and parts[0] == '3':
                try:
                    data_eq.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_eq.append(row)


    # Adjust column names to the file
    columns = ["RUN" ,"PT", "TY", "LAB", "I", "L2NORM", "v", "z", "stability"]
    df_eq = pd.DataFrame(data_eq, columns=columns)

    hopf = df_eq[df_eq["TY"] == 3]
    saddle_nodes = df_eq[df_eq["TY"] == 2]

    original_eq = df_eq[df_eq["RUN"] != 2]
    stable_eq_og = original_eq[original_eq["stability"] == 2]
    stable_eq_og = stable_eq_og.copy()
    stable_eq_og.loc[stable_eq_og['TY'] == 2, ["I","v"]] = np.nan
    stable_eq_og.loc[stable_eq_og['TY'] == 3, ["I","v"]] = np.nan

    unstable_eq_og = original_eq[original_eq['stability'] != 2]
    unstable_eq_og = unstable_eq_og.copy()
    unstable_eq_og.loc[unstable_eq_og['TY'] == 2, ["I","v"]] = np.nan
    unstable_eq_og.loc[unstable_eq_og['TY'] == 3, ["I","v"]] = np.nan

    fs = df_eq[df_eq["RUN"] == 2]
    stable_eq_fs = fs[fs["stability"] == 2]
    stable_eq_fs = stable_eq_fs.copy()
    stable_eq_fs.loc[stable_eq_fs['TY'] == 2, ["I","v"]] = np.nan
    stable_eq_fs.loc[stable_eq_fs['TY'] == 3, ["I","v"]] = np.nan

    unstable_eq_fs = fs[fs['stability'] != 2]
    unstable_eq_fs = unstable_eq_fs.copy()
    unstable_eq_fs.loc[unstable_eq_fs['TY'] == 2, ["I","v"]] = np.nan
    unstable_eq_fs.loc[unstable_eq_fs['TY'] == 3, ["I","v"]] = np.nan

    plt.figure(figsize=(7,5))
    plt.plot(stable_eq_og["I"], stable_eq_og["v"], label="stable equilibria", color="blue")
    plt.plot(unstable_eq_og["I"], unstable_eq_og["v"], label="unstable equilibria", color="blue", ls='--')
    plt.plot(stable_eq_fs["I"], stable_eq_fs["v"], label="stable FS", color="red")
    plt.plot(unstable_eq_fs["I"], unstable_eq_fs["v"], label="unstable FS", color="red", ls='--')
    
    #plot the hopf point
    if not hopf.empty:
        plt.scatter(hopf["I"], hopf["v"], color="red", marker="o", s=20, label="Hopf Bifurcation", zorder=10)
    if not saddle_nodes.empty:
        plt.scatter(saddle_nodes["I"], saddle_nodes["v"], color="black", marker="x", s=20, label="Saddle Node Bifurcation", zorder=10)

    plt.xlabel("I (non-dimensional stimulus current)")
    plt.ylabel("v (non-dimensional voltage)")
    plt.title(figurename)
    plt.legend()
    plt.grid(True)
    plt.savefig(figurename)
    plt.close()


# todo: I just want to remove the first/last point from each branch when it switches from stable to unstable and otherwise.
# tonights job...