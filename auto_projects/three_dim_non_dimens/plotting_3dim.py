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
            if parts and parts[0] == '-4':
                try:
                    data_periodic.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_periodic.append(row)
            if parts and parts[0] == '-2':
                try:
                    data_periodic.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_periodic.append(row)


    # Adjust column names to your file
    columns = ["RUN" ,"PT", "TY", "LAB", "I", "L2NORM", "v", "y", "z", "stability"]
    periodic = ["RUN" ,"PT", "TY", "LAB", "I", "L2NORM", "max v", "max y", "max z", "period", "min v", "stability"]
    df_eq = pd.DataFrame(data_eq, columns=columns)
    df_periodic = pd.DataFrame(data_periodic, columns=periodic)
    df_periodic.iloc[:4] = np.nan

    hopf = df_eq[df_eq["TY"] == 3]
    saddle_nodes = df_eq[df_eq["TY"] == 2]

    stable_eq = df_eq[df_eq["stability"] == 3]
    stable_eq = stable_eq.copy()
    stable_eq.loc[stable_eq['TY'] == 2, ["I","v"]] = np.nan

    unstable_eq = df_eq[df_eq['stability'] != 3]
    unstable_eq = unstable_eq.copy()
    unstable_eq.loc[unstable_eq['TY'] == 2, ["I","v"]] = np.nan
    
    # I want to edit this so that if there is a difference between the start and stuff
    # rather than looking for the branch point because it is a little off, maybe we make like after the branch point the next five are nan or somethig??
    unstable_limit_cycle = df_periodic[df_periodic['stability'] != 3]
    if filename == "b.one_par_type1":
        unstable_limit_cycle = unstable_limit_cycle.copy()
        unstable_limit_cycle.loc[unstable_limit_cycle['PT'] == 1, ["I", "max v", "min v"]] = np.nan

    stable_limit_cycle = df_periodic[df_periodic['stability'] == 3]
    if filename == "b.one_par_type2":
        stable_limit_cycle = stable_limit_cycle.copy()

        idxs = stable_limit_cycle.index[stable_limit_cycle['TY'] == 5].tolist()

        for idx in idxs:
            pos = stable_limit_cycle.index.get_loc(idx)
            affected = stable_limit_cycle.iloc[pos : pos+20].index

            stable_limit_cycle.loc[affected, ["I", "max v", "min v"]] = np.nan

    # Plot min v and max v vs I
    plt.figure(figsize=(7,5))
    plt.plot(stable_eq["I"], stable_eq["v"], label="stable equilibria", color="blue")
    plt.plot(unstable_eq["I"], unstable_eq["v"], label="unstable equilibria", color="red", ls = '--')
    plt.plot(stable_limit_cycle["I"], stable_limit_cycle["min v"], color='black', label="stable periodic branch")
    plt.plot(stable_limit_cycle["I"], stable_limit_cycle['max v'], color='black')
    plt.plot(unstable_limit_cycle["I"], unstable_limit_cycle["min v"], color='black', label="unstable periodic branch", ls = '--')
    plt.plot(unstable_limit_cycle["I"], unstable_limit_cycle['max v'], color='black', ls = '--')

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


def plot_two_param(filename, figurename, axes, ylabel):
    data_saddle = []
    data_hopf = []
    with open(filename) as f:
        for line in f:
            parts = line.strip().split()
            if parts and (filename == "b.two_par_cusp") and parts[0] == '2':
                try:
                    data_saddle.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_saddle.append(row)
            if parts and parts[0] == '1':
                try:
                    data_hopf.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_hopf.append(row)
            if parts and (filename == "b.two_par_hopf" or filename == "b.two_par_hopf_alt") and parts[0] == '2':
                try:
                    data_hopf.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_hopf.append(row)
            if parts and (filename == "b.two_par_low_eps") and parts[0]== '2':
                try:
                    data_hopf.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_hopf.append(row)


    # Adjust column names to your file
    columns = ["RUN" ,"PT", "TY", "LAB", "I", "L2NORM", "v", "y", "z", "gsub", 'stability']
    periodic = ["RUN" ,"PT", "TY", "LAB", "I", "L2NORM", "v", "y", "z", "gsub", "stability"]
    df_saddle = pd.DataFrame(data_saddle, columns=columns)
    df_hopf = pd.DataFrame(data_hopf, columns=periodic)
    
    bt = pd.DataFrame()
    if not df_saddle.empty:
        bt = df_saddle[(df_saddle["TY"] == -21) | (df_saddle["TY"] == -31)]
    if not df_hopf.empty:
        bt = df_hopf[(df_hopf["TY"] == -21) | (df_hopf["TY"] == -31)]
    cusp = df_saddle[df_saddle["TY"] == -22]
    gh = df_hopf[df_hopf["TY"] == -32]

    # Plot min v and max v vs I
    plt.figure(figsize=(7,5))
    if not df_saddle.empty:
        plt.plot(df_saddle["I"], df_saddle["gsub"], label="curve of saddle nodes", color="red")
    if not bt.empty:
        plt.scatter(bt["I"], bt["gsub"], color="black", label="BT", zorder=10)
    if not cusp.empty:
        plt.scatter(cusp["I"], cusp["gsub"], color='purple', label="Cusp Point", zorder=10)
    if not df_hopf.empty:
        df_hopf.loc[df_hopf['TY'] == -34, ["I","v"]] = np.nan
        plt.plot(df_hopf["I"], df_hopf["gsub"], color='blue', label="curve of hopf bifurcations")
    if not gh.empty:
        plt.scatter(gh["I"], gh["gsub"], color="orange", label="GH", zorder=10)
    plt.xlabel("I (non-dimensional stimulus current)")
    plt.ylabel(ylabel)
    plt.title(figurename)
    plt.legend()
    plt.grid(True)
    # to change for 'zoomed in' plots etc.
    plt.ylim(axes[2], axes[3])
    plt.xlim(axes[0], axes[1])
    plt.savefig(figurename)
    plt.close()

def plot_one_param_zoom(filename, figurename, zoom_region, zoom_periodic=False):
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
            if parts and parts[0] == '-4':
                try:
                    data_periodic.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_periodic.append(row)
            if parts and parts[0] == '-2':
                try:
                    data_periodic.append([float(x) for x in parts])
                except ValueError:
                    row = []
                    for x in parts:
                        try:
                            row.append(float(x))
                        except ValueError:
                            row.append(x)
                    data_periodic.append(row)


    # Adjust column names to your file
    columns = ["RUN" ,"PT", "TY", "LAB", "I", "L2NORM", "v", "y", "z", "stability"]
    periodic = ["RUN" ,"PT", "TY", "LAB", "I", "L2NORM", "max v", "max y", "max z", "period", "min v", "stability"]
    df_eq = pd.DataFrame(data_eq, columns=columns)
    df_eq = df_eq.copy()
    df_eq.loc[(df_eq['I'] < zoom_region[0]) | (df_eq['I'] > zoom_region[1]), ["I","v"]] = np.nan
    df_periodic = pd.DataFrame(data_periodic, columns=periodic)
    df_periodic = df_periodic.copy()
    df_periodic.loc[(df_periodic['I'] < zoom_region[0]) | (df_periodic['I'] > zoom_region[1]), ["I","max v", "min v"]] = np.nan
    if zoom_periodic > 2:
        df_eq = df_eq.copy()
        df_eq.loc[(df_eq['I'] < zoom_region[0]) | (df_eq['I'] > zoom_region[1]), ["I","v"]] = np.nan
        df_periodic = df_periodic.copy()
        df_periodic.loc[(df_periodic['min v'] < zoom_region[2]) | (df_periodic['max v'] > zoom_region[3]), ["I","max v", "min v"]] = np.nan

    hopf = df_eq[df_eq["TY"] == 3]
    saddle_nodes = df_eq[df_eq["TY"] == 2]

    stable_eq = df_eq[df_eq["stability"] == 3]
    stable_eq = stable_eq.copy()
    stable_eq.loc[stable_eq['TY'] == 2, ["I","v"]] = np.nan

    unstable_eq = df_eq[df_eq['stability'] != 3]
    unstable_eq = unstable_eq.copy()
    unstable_eq.loc[unstable_eq['TY'] == 2, ["I","v"]] = np.nan

    # now the first point is labelled unstable, but I don't think it is meant to
    unstable_limit_cycle = df_periodic[df_periodic['stability'] != 3]
    if filename == "b.one_par_type1":
        unstable_limit_cycle = unstable_limit_cycle.copy()
        unstable_limit_cycle.loc[unstable_limit_cycle['PT'] == 1, ["I", "max v", "min v"]] = np.nan
    
    stable_limit_cycle = df_periodic[df_periodic['stability'] == 3]

    # Plot min v and max v vs I
    plt.figure(figsize=(7,5))
    plt.plot(stable_eq["I"], stable_eq["v"], label="stable equilibria", color="blue")
    plt.plot(unstable_eq["I"], unstable_eq["v"], label="unstable equilibria", color="red", ls = '--')
    plt.plot(stable_limit_cycle["I"], stable_limit_cycle["min v"], color='black')
    plt.plot(stable_limit_cycle["I"], stable_limit_cycle['max v'], color='black', label="stable periodic branch")
    plt.plot(unstable_limit_cycle["I"], unstable_limit_cycle["min v"], color='black', label="unstable periodic branch", ls = '--')
    plt.plot(unstable_limit_cycle["I"], unstable_limit_cycle['max v'], color='black', ls = '--')

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

# if __name__ == "__main__":