import matplotlib.pyplot as plt
from analyses import *
import copy 


def plot_lambda_hist(PATH_ANALYSIS, xvg_data, coord2lambda_dict, lambda_ref, rows=20, cols=5, quality='Debug'):
    """Plot histogram of lambda values"""
    # Set up the grid size
    # Create a figure with subplots
    # Adjusted figure size for 5x20 layout
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))
    #print(xvg_data.coordids)
    coordid = xvg_data.coordids[0]
    #print(coordid)
    index_coordid = copy.deepcopy(xvg_data.coordids[0])  # Create a copy of the coordid list to iterate through
    #coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    if quality == 'Debug':
        fig.suptitle(PATH_ANALYSIS, fontsize=16, y=0.995)  # Moved suptitle up
    elif quality == 'Publication':
        fig.suptitle("Lambda distributions", fontsize=16,
                     y=0.995)  # Moved suptitle up

    # Adjust layout to leave space for suptitle
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)

    for i in range(rows):
        for j in range(cols):
            index = coordid - index_coordid
            #print(index)

            if index < len(xvg_data.coordids):
                ax = axes[i, j]
                data = xvg_data[coordid]
                ax.hist(data[:, 1], bins=500)

                ax.set_xlim(-0.125, 1.125)
                # Decide whether to use residue name or coordid in title
                if quality == 'Debug':
                    ax.set_title(
                        f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')
                elif quality == 'Publication':
                    ax.set_title(
                        f'{lambda_ref.iloc[coordid-1]["resname"]}_{lambda_ref.iloc[coordid-1]["resid"]}')
            else:
                continue

            coordid += 1

    for ax in axes.flat[coordid-1:]:
        ax.remove()

# Show the plot
    return plt


def plot_protonation_timeseries(PATH_ANALYSIS, time, xvg_data, coord2lambda_dict, lambda_ref, rows=20, cols=5, quality='Debug', npz_output=False):
    """"Plot protonation time-series"""
    # Set up the grid size
    # Create a figure with subplots
    # Adjusted figure size for 5x20 layout
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))

    coordid = xvg_data.coordids[0]
    #print(coordid)
    index_coordid = copy.deepcopy(xvg_data.coordids[0])
    #coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    if quality == 'Debug':
        fig.suptitle(PATH_ANALYSIS, fontsize=16, y=0.995)  # Moved suptitle up
    elif quality == 'Publication':
        fig.suptitle("Protonation time-series", fontsize=16,
                     y=0.995)  # Moved suptitle up

    # Adjust layout to leave space for suptitle
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)

    for i in range(rows):
        for j in range(cols):
            index = coordid - index_coordid

            if index < len(xvg_data.coordids):
                ax = axes[i, j]

                if lambda_ref.iloc[index]['resname'] == "HSPT":  # If histidine
                    # If this is histidine lambda1
                    if lambda_ref.iloc[index]['resid'] != prv_resid:

                        coordids = [coordid, coordid+1, coordid+2]

                        res_prot_ts = get_histidine_protonation_timeseries(
                            coordids, xvg_data)

                        prv_resid = lambda_ref.iloc[index]['resid']
                    else:
                        res_prot_ts = get_histidine_protonation_timeseries(
                            coordids, xvg_data)

                else:
                    res_prot_ts = get_protonation_timeseries(
                        coordid, xvg_data)

                if npz_output:
                    # Save output in npz file named after the residue
                    resname = lambda_ref.iloc[index]['resname']
                    resid = lambda_ref.iloc[index]['resid']
                    np.savez(
                        f"{PATH_ANALYSIS}/{resname}_{resid}_protonation_timeseries.npz", res_prot_ts=res_prot_ts)

                ax.plot(res_prot_ts, label="MD1")
                ax.set_ylim(-0.1, 1.1)

                # Set xticks and labels aaccording to the simulation time
                length = round(time/10000)*10000
                xticks = np.arange(0, length, length/3)
                xticks = np.round(xticks/10000)*10000
                xticks = np.concatenate((xticks, [length]))
                ax.set_xticks(xticks.astype(int))
                ax.set_xticklabels((xticks/1000).astype(int))

                if quality == 'Debug':
                    ax.set_title(
                        f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')
                elif quality == 'Publication':
                    ax.set_title(
                        f'{lambda_ref.iloc[coordid-1]["resname"]}_{lambda_ref.iloc[coordid-1]["resid"]}')

            else:
                continue

            coordid += 1

    for ax in axes.flat[coordid-1:]:
        ax.remove()

    return plt


def plot_protonation_convergence(PATH_ANALYSIS, time, xvg_data_list: List[XVGData], coord2lambda_dict, lambda_ref, rows=20, cols=5, quality='Debug'):
    """Plot protonation avg and standard error time-series. Just two replicas supported - add get_protfrac_ts if more than 2"""
    # Set up the grid size
    # Create a figure with subplots
    # Adjusted figure size for 5x20 layout
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    if quality == 'Debug':
        fig.suptitle(PATH_ANALYSIS, fontsize=16, y=0.995)  # Moved suptitle up
    elif quality == 'Publication':
        fig.suptitle("Protonation convergence", fontsize=16,
                     y=0.995)  # Moved suptitle up

    # Adjust layout to leave space for suptitle
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)

    for i in range(rows):
        for j in range(cols):
            list_residues = []

            index = coordid-1

            if index < lambda_ref.shape[0]:
                ax = axes[i, j]

                if lambda_ref.iloc[index]['resname'] == "HSPT":  # If histidine
                    # If this is histidine lambda1
                    if lambda_ref.iloc[index]['resid'] != prv_resid:

                        coordids = [coordid, coordid+1, coordid+2]

                        for xvg_data in xvg_data_list:
                            list_residues.append(get_histidine_protonation_timeseries(
                                coordids, xvg_data))

                        prv_resid = lambda_ref.iloc[index]['resid']
                    else:
                        for xvg_data in xvg_data_list:
                            list_residues.append(get_histidine_protonation_timeseries(
                                coordids, xvg_data))
                else:
                    for xvg_data in xvg_data_list:
                        list_residues.append(get_protonation_timeseries(
                            coordid, xvg_data))

                # min_length = min(len(res_prot_ts1), len(res_prot_ts2))
                min_length = min(map(len, list_residues))
                total_protarray = np.vstack(
                    [res_array[:min_length] for res_array in list_residues])
                total_protonse = np.std(
                    total_protarray, axis=0) / np.sqrt(len(total_protarray))
                ax.plot(np.mean(total_protarray, axis=0))
                ax.fill_between(np.arange(len(total_protarray[0, :])), np.mean(
                    total_protarray, axis=0) - total_protonse, np.mean(total_protarray, axis=0) + total_protonse, alpha=0.5)

                # Set xticks and labels aaccording to the simulation time
                length = round(time/10000)*10000
                xticks = np.arange(0, length, length/3)
                xticks = np.round(xticks/10000)*10000
                xticks = np.concatenate((xticks, [length]))
                ax.set_xticks(xticks.astype(int))
                ax.set_xticklabels((xticks/1000).astype(int))

                ax.set_ylim(-0.1, 1.1)

                # Decide whether to use residue name or coordid in title
                if quality == 'Debug':
                    ax.set_title(
                        f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')
                elif quality == 'Publication':
                    ax.set_title(
                        f'{lambda_ref.iloc[coordid-1]["resname"]}_{lambda_ref.iloc[coordid-1]["resid"]}')

            else:
                continue

            coordid += 1

    for ax in axes.flat[coordid-1:]:
        ax.remove()

    return plt


def plot_protonation_fraction(xvg_data_list: List[XVGData], lambda_ref, rows=20, cols=5, npz_output=False):
    """Plot protonation fractions avg and se"""
    # Set up the grid size
    # Create a figure with subplots
    # Adjusted figure size for 5x20 layout
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))
    fig.tight_layout(pad=2.0)  # Optional: adjust padding between plots

    coordid = 1
    prv_resid = 0

    # Create folder for npz protonation fraction files
    if npz_output:
        npz_dir = f"npz_protfrac"
        if not os.path.exists(npz_dir):
            os.makedirs(npz_dir)

    # Generate and plot data for each subplot
    for i in range(rows):
        for j in range(cols):
            index = coordid-1

            if index < lambda_ref.shape[0]:
                ax = axes[i, j]

                # Generate some data; here, we're using sine waves with different frequencies

                if lambda_ref.iloc[index]['resname'] == "HSPT":  # If histidine
                    # If this is histidine lambda1
                    if lambda_ref.iloc[index]['resid'] != prv_resid:

                        coordids = [coordid, coordid+1, coordid+2]

                        proton_avg, proton_se = get_histidine_statistics(
                            coordids, xvg_data_list)

                        prv_resid = lambda_ref.iloc[index]['resid']
                    else:
                        proton_avg, proton_se = get_histidine_statistics(
                            coordids, xvg_data_list)

                else:
                    proton_avg, deproton_avg, proton_se, deproton_se = get_statistics(
                        coordid, xvg_data_list)

                bars = ax.bar(['Protonated'], [proton_avg], yerr=proton_se,
                              capsize=5, linewidth=1, edgecolor='black')
                ax.bar_label(
                    bars, labels=[f"{proton_avg:.2f} ± {proton_se:.2f}"])

                ax.set_ylim(0, 1.2)
                ax.set_title(
                    f'{lambda_ref.iloc[coordid-1]["resname"]}_{lambda_ref.iloc[coordid-1]["resid"]}')

                ax.set_yticks([0, 0.5, 1])
                ax.set_xticks([])

                # Save output in npz file named after the residue
                if npz_output:
                    resname = lambda_ref.iloc[index]['resname']
                    resid = lambda_ref.iloc[index]['resid']
                    np.savez(f"{npz_dir}/{resname}_{resid}_protonation_fraction.npz",
                             res_prot_avg=proton_avg, res_prot_se=proton_se)

            else:
                continue

            coordid += 1

    # Show the plot
    for ax in axes.flat[coordid-1:]:
        ax.remove()

    return plt


def single_residue_convergence(coordid, xvg_data_list: List[XVGData], lambda_ref, title="Constant-pH MD"):
    """THIS WORKS ONLY FOR NON HISTIDINES! Plot convergence of single residue. Just two replicas supported - add res_fracX if more."""

    res_fractions = []

    for xvg_data in xvg_data_list:
        res_fractions.append(get_protonation_timeseries(coordid, xvg_data))
    # res_frac2 = get_protonation_timeseries(coordid, PATH_ANALYSIS[1], xvg_data)

    min_length = min(map(len, res_fractions))

    total_protarray = np.vstack(
        [res_frac[:min_length] for res_frac in res_fractions])

    total_protonse = np.std(total_protarray, axis=0) / \
        np.sqrt(len(total_protarray))

    plt.plot(np.mean(total_protarray, axis=0),
             label=f'{lambda_ref.iloc[coordid-1]["resname"]}_{lambda_ref.iloc[coordid-1]["resid"]}')
    plt.fill_between(np.arange(len(total_protarray[0, :])), np.mean(
        total_protarray, axis=0) - total_protonse, np.mean(total_protarray, axis=0) + total_protonse, alpha=0.5)

    plt.ylim(0, 1)
    plt.ylabel("Protonation fraction")
    plt.xlabel("Time (ps)")
    plt.legend()
    plt.title(title)
    return plt
