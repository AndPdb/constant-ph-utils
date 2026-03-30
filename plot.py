import matplotlib.pyplot as plt
from analyses import *
import copy
from collections import OrderedDict


def _group_coordids_by_resname(coordids, lambda_ref):
    """Group coordids by their three-letter residue name.

    Returns an OrderedDict: {resname: [coordid, ...]} preserving
    the order in which each residue type first appears.
    """
    groups = OrderedDict()
    for cid in coordids:
        resname = lambda_ref.iloc[cid - 1]["resname"]
        groups.setdefault(resname, []).append(cid)
    return groups


def plot_lambda_hist(xvg_data, coord2lambda_dict, lambda_ref, rows=20, cols=5, quality='Debug'):
    """Plot histogram of lambda values"""
    # Set up the grid size
    # Create a figure with subplots
    # Adjusted figure size for 5x20 layout
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))
    # print(xvg_data.coordids)
    coordid = xvg_data.coordids[0]
    # print(coordid)
    # Create a copy of the coordid list to iterate through
    index_coordid = copy.deepcopy(xvg_data.coordids[0])
    # coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    if quality == 'Debug':
        fig.suptitle(xvg_data.analysis_dir, fontsize=16,
                     y=0.995)  # Moved suptitle up
    elif quality == 'Publication':
        fig.suptitle("Lambda distributions", fontsize=16,
                     y=0.995)  # Moved suptitle up

    # Adjust layout to leave space for suptitle
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)

    for i in range(rows):
        for j in range(cols):
            index = coordid - index_coordid

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


def plot_protonation_timeseries(time, xvg_data, coord2lambda_dict, lambda_ref, rows=20, cols=5, quality='Debug', npz_output=False):
    """"Plot protonation time-series"""
    # Set up the grid size
    # Create a figure with subplots
    # Adjusted figure size for 5x20 layout
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))

    coordid = xvg_data.coordids[0]
    # print(coordid)
    index_coordid = copy.deepcopy(xvg_data.coordids[0])
    # coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    if quality == 'Debug':
        fig.suptitle(xvg_data.analysis_dir, fontsize=16,
                     y=0.995)  # Moved suptitle up
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
                        f"{xvg_data.analysis_dir}/{resname}_{resid}_protonation_timeseries.npz", res_prot_ts=res_prot_ts)

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


def plot_protonation_convergence(PATH_ANALYSIS, time, xvg_data_list: List[XVGData], coord2lambda_dict, lambda_ref, chain_mapping={}, rows=20, cols=5, quality='Debug'):
    """Plot protonation avg and standard error time-series. Just two replicas supported - add get_protfrac_ts if more than 2"""
    # Set up the grid size
    # Create a figure with subplots
    # Adjusted figure size for 5x20 layout
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))
    dimensions = len(xvg_data_list[0].coordids)
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
    length = round(time/10000)*10000
    xticks = np.arange(0, length, length/3)
    xticks = np.round(xticks/10000)*10000
    xticks = np.concatenate((xticks, [length]))

    for i in range(rows):
        for j in range(cols):
            list_residues = []

            index = coordid-1

            if index < dimensions:
                ax = axes[i, j]

                if lambda_ref.iloc[index]['resname'] == "HSPT":  # If histidine
                    # If this is histidine lambda1
                    if lambda_ref.iloc[index]['resid'] != prv_resid:

                        coordids = [coordid, coordid+1, coordid+2]

                        for xvg_data in xvg_data_list:
                            list_residues.append(get_histidine_protonation_timeseries(
                                coordids, xvg_data, chain_mapping))

                        prv_resid = lambda_ref.iloc[index]['resid']
                    else:
                        for xvg_data in xvg_data_list:
                            list_residues.append(get_histidine_protonation_timeseries(
                                coordids, xvg_data, chain_mapping))
                else:
                    for xvg_data in xvg_data_list:
                        list_residues.append(get_protonation_timeseries(
                            coordid, xvg_data, chain_mapping))

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


def plot_protonation_fraction(xvg_data_list: List[XVGData], lambda_ref,
                              chain_mapping={}, npz_output=False):
    """Plot protonation fractions, one figure per residue type.
    All residues of the same type are shown as bars on a single axis.
    """
    all_coordids = list(range(1, len(xvg_data_list[0].coordids) + 1))
    groups = _group_coordids_by_resname(all_coordids, lambda_ref)
    figures = {}

    if npz_output:
        npz_dir = "npz_protfrac"
        if not os.path.exists(npz_dir):
            os.makedirs(npz_dir)

    for resname, cids in groups.items():
        labels = []
        avgs = []
        ses = []
        prv_resid = 0
        seen_resids = set()

        for coordid in cids:
            index = coordid - 1
            resid = lambda_ref.iloc[index]['resid']

            if lambda_ref.iloc[index]['resname'] == "HSPT":
                if resid in seen_resids:
                    continue
                seen_resids.add(resid)
                if resid != prv_resid:
                    hist_coordids = [coordid, coordid + 1, coordid + 2]
                    prv_resid = resid
                proton_avg, proton_se = get_histidine_statistics(
                    hist_coordids, xvg_data_list, chain_mapping)
            else:
                proton_avg, deproton_avg, proton_se, deproton_se = \
                    get_statistics(coordid, xvg_data_list, chain_mapping)

            labels.append(f"{resname}_{resid}")
            avgs.append(proton_avg)
            ses.append(proton_se)

            if npz_output:
                rn = lambda_ref.iloc[index]['resname']
                np.savez(f"{npz_dir}/{rn}_{resid}_protonation_fraction.npz",
                         res_prot_avg=proton_avg, res_prot_se=proton_se)

        n = len(avgs)
        bar_color = "#4472C4"
        fig, ax = plt.subplots(figsize=(max(3, n * 1.2), 4))
        x = np.arange(n)
        bars = ax.bar(x, avgs, yerr=ses, capsize=5, width=0.5,
                      color=bar_color, edgecolor=bar_color,
                      error_kw=dict(lw=1.2, capthick=1.2))
        ax.bar_label(bars,
                     labels=[f"{a:.2f} ± {s:.2f}" for a, s in zip(avgs, ses)],
                     padding=5, fontsize=8)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=10)
        ax.set_ylabel("Protonation Fraction", fontsize=11)
        ax.set_ylim(0.0, 1.15)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', which='both', labelsize=10)
        fig.tight_layout()

        figures[resname] = fig

    return figures


def single_residue_convergence(coordid, xvg_data_list: List[XVGData], lambda_ref, chain_mapping={}, title="Constant-pH MD"):
    """THIS WORKS ONLY FOR NON HISTIDINES! Plot convergence of single residue. Just two replicas supported - add res_fracX if more."""

    res_fractions = []

    for xvg_data in xvg_data_list:
        res_fractions.append(get_protonation_timeseries(
            coordid, xvg_data, chain_mapping=chain_mapping))

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
