import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_one_param(filename, figurenames, zoom_type2=False, zoom_type3=False):
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
    stable_eq_og.loc[stable_eq_og['TY'].isin([2, 3, 6, 9, -4]), ["I","v"]] = np.nan
    # stable_eq_og.loc[stable_eq_og['TY'] == 3, ["I","v"]] = np.nan

    unstable_eq_og = original_eq[original_eq['stability'] != 2]
    unstable_eq_og = unstable_eq_og.copy()
    unstable_eq_og.loc[unstable_eq_og['TY'].isin([2, 3, 6, 9, -4]), ["I","v"]] = np.nan
    # unstable_eq_og.loc[unstable_eq_og['TY'] == 3, ["I","v"]] = np.nan

    fs = df_eq[df_eq["RUN"] == 2]
    stable_eq_fs = fs[fs["stability"] == 2]
    stable_eq_fs = stable_eq_fs.copy()
    stable_eq_fs.loc[stable_eq_fs['TY'].isin([2, 3, 6, 9, -4]), ["I","v"]] = np.nan
    # stable_eq_fs.loc[stable_eq_fs['TY'] == 3, ["I","v"]] = np.nan

    unstable_eq_fs = fs[fs['stability'] != 2]
    unstable_eq_fs = unstable_eq_fs.copy()
    unstable_eq_fs.loc[unstable_eq_fs['TY'].isin([2, 3, 6, 9, -4]), ["I","v"]] = np.nan
    # unstable_eq_fs.loc[unstable_eq_fs['TY'] == 3, ["I","v"]] = np.nan

    plt.figure(figsize=(7,5))
    plt.plot(stable_eq_og["I"]/(100*20), stable_eq_og["v"]/100, label="stable equilibria", color="blue")
    plt.plot(unstable_eq_og["I"]/(100*20), unstable_eq_og["v"]/100, label="unstable equilibria", color="blue", ls='--')
    plt.plot(stable_eq_fs["I"]/(100*20), stable_eq_fs["v"]/100, label="stable FS", color="red")
    plt.plot(unstable_eq_fs["I"]/(100*20), unstable_eq_fs["v"]/100, label="unstable FS", color="red", ls='--')
    
    #plot the hopf point
    if not hopf.empty:
        plt.scatter(hopf["I"]/(100*20), hopf["v"]/100, color="red", marker="o", s=20, label="Hopf Bifurcation", zorder=10)
    if not saddle_nodes.empty:
        plt.scatter(saddle_nodes["I"]/(100*20), saddle_nodes["v"]/100, color="black", marker="x", s=20, label="Saddle Node Bifurcation", zorder=10)

    plt.xlabel("I (non-dimensional stimulus current)")
    plt.ylabel("v (non-dimensional voltage)")
    plt.title(figurenames[0])
    plt.legend()
    plt.grid(True)
    plt.savefig(figurenames[0])

    # now, we also want to plot in z.
    plt.figure(figsize=(7,5))
    plt.plot(stable_eq_og["I"]/(100*20), stable_eq_og["z"], label="stable equilibria", color="blue")
    plt.plot(unstable_eq_og["I"]/(100*20), unstable_eq_og["z"], label="unstable equilibria", color="blue", ls='--')
    plt.plot(stable_eq_fs["I"]/(100*20), stable_eq_fs["z"], label="stable FS", color="red")
    plt.plot(unstable_eq_fs["I"]/(100*20), unstable_eq_fs["z"], label="unstable FS", color="red", ls='--')
    
    #plot the hopf point
    if not hopf.empty:
        plt.scatter(hopf["I"]/(100*20), hopf["z"], color="red", marker="o", s=20, label="Hopf Bifurcation", zorder=10)
    if not saddle_nodes.empty:
        plt.scatter(saddle_nodes["I"]/(100*20), saddle_nodes["z"], color="black", marker="x", s=20, label="Saddle Node Bifurcation", zorder=10)

    plt.xlabel("I (non-dimensional stimulus current)")
    plt.ylabel("z (non-dimensional gating term)")
    plt.title(figurenames[1])
    plt.legend()
    plt.grid(True)
    plt.xlim(0, 0.08)
    plt.ylim(-1, 1)
    plt.savefig(figurenames[1])
    plt.close()

    if zoom_type2:
        plt.figure(figsize=(7,5))
        plt.plot(stable_eq_og["I"]/(100*20), stable_eq_og["z"], label="stable equilibria", color="blue")
        plt.plot(unstable_eq_og["I"]/(100*20), unstable_eq_og["z"], label="unstable equilibria", color="blue", ls='--')
        plt.plot(stable_eq_fs["I"]/(100*20), stable_eq_fs["z"], label="stable FS", color="red")
        plt.plot(unstable_eq_fs["I"]/(100*20), unstable_eq_fs["z"], label="unstable FS", color="red", ls='--')
        plt.plot([45/(100*20), 45/(100*20)], [-0.3,1], '--', color='black', label="One Transient Spike Observed")

        # colour the region
        unstable_eq_fs['between_line'] = unstable_eq_fs['z'].clip(lower=0)

        I1 = unstable_eq_fs['I'].to_numpy()
        z1 = unstable_eq_fs['between_line'].to_numpy()

        order1 = np.argsort(I1)
        I1_sorted = I1[order1]
        z1_sorted = z1[order1]

        # Stable curve
        I2 = stable_eq_og['I'].to_numpy()
        z2 = stable_eq_og['z'].to_numpy()

        order2 = np.argsort(I2)
        I2_sorted = I2[order2]
        z2_sorted = z2[order2]

        # Combined x-axis in I-space
        xfill = np.sort(np.concatenate([I1_sorted, I2_sorted]))

        # Interpolate z onto this I grid
        y1fill = np.interp(xfill, I1_sorted, z1_sorted)
        y2fill = np.interp(xfill, I2_sorted, z2_sorted)
        
        plt.fill_between(xfill, y1fill, y2fill, where=y1fill > y2fill, interpolate=True, color='crimson', alpha=0.2)
        plt.fill_between(xfill, y1fill, y2fill, where=y1fill < y2fill, interpolate=True, color='blue', alpha=0.2)

        #plt.fill_between(stable_eq_og["I"]/(100*20), stable_eq_og["z"], 0, alpha=0.2, color="grey")
        
        #plot the hopf point
        if not hopf.empty:
            plt.scatter(hopf["I"]/(100*20), hopf["z"], color="red", marker="o", s=20, label="Hopf Bifurcation", zorder=10)
        if not saddle_nodes.empty:
            plt.scatter(saddle_nodes["I"]/(100*20), saddle_nodes["z"], color="black", marker="x", s=20, label="Saddle Node Bifurcation", zorder=10)

        plt.xlabel("I (non-dimensional stimulus current)")
        plt.ylabel("z (non-dimensional gating term)")
        plt.title(figurenames[2])
        plt.legend()
        plt.grid(True)
        plt.xlim(0.0, 0.025)
        plt.ylim(-0.3, 1)
        plt.savefig(figurenames[2])
        plt.close()
    
    if zoom_type3: 
        plt.figure(figsize=(7,5))
        plt.plot(stable_eq_og["I"]/(100*20), stable_eq_og["z"], label="stable equilibria", color="blue")
        plt.plot(unstable_eq_og["I"]/(100*20), unstable_eq_og["z"], label="unstable equilibria", color="blue", ls='--')
        plt.plot(stable_eq_fs["I"]/(100*20), stable_eq_fs["z"], label="stable FS", color="red")
        plt.plot(unstable_eq_fs["I"]/(100*20), unstable_eq_fs["z"], label="unstable FS", color="red", ls='--')
        plt.plot([80/(100*20), 80/(100*20)], [0.05,0.2], '--', color='black', label="Two Transient Spikes Observed")


        #plot the hopf point
        if not hopf.empty:
            plt.scatter(hopf["I"]/(100*20), hopf["z"], color="red", marker="o", s=20, label="Hopf Bifurcation", zorder=10)
        if not saddle_nodes.empty:
            plt.scatter(saddle_nodes["I"]/(100*20), saddle_nodes["z"], color="black", marker="x", s=20, label="Saddle Node Bifurcation", zorder=10)

        plt.xlabel("I (non-dimensional stimulus current)")
        plt.ylabel("z (non-dimensional gating term)")
        plt.title(figurenames[2])
        plt.legend()
        plt.grid(True)
        plt.xlim(0.035, 0.045)
        plt.ylim(0.05, 0.2)
        plt.savefig(figurenames[2])
        plt.close()
# todo: I just want to remove the first/last point from each branch when it switches from stable to unstable and otherwise.
# tonights job...