import matplotlib.pyplot as plt
from analyses import *

def plot_lambda_hist(PATH_ANALYSIS, xvg_data, coord2lambda_dict, lambda_ref, rows=20, cols=5):
    """Plot histogram of lambda values"""
        # Set up the grid size
    #rows = 20  # 10x10 grid for 100 plots
    #cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    fig.suptitle(PATH_ANALYSIS, fontsize=16, y=0.995)  # Moved suptitle up
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)  # Adjust layout to leave space for suptitle

    for i in range(rows):
        for j in range(cols):
            index = coordid-1
            #print(f"Coordid {coordid}")

            if index < lambda_ref.shape[0]:
                ax = axes[i, j]
                data = xvg_data.get_coord_xvg(coordid, PATH_ANALYSIS)
                ax.hist(data[:,1], bins=500);

                ax.set_xlim(-0.125,1.125)
                ax.set_title(f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')

            else:
                continue
            
            coordid += 1

    for ax in axes.flat[coordid-1:]:
        ax.remove()

# Show the plot
    return plt

def plot_protonation_timeseries(PATH_ANALYSIS, time, xvg_data, coord2lambda_dict, lambda_ref, rows=20, cols=5, npz_output=False):
    """"Plot protonation time-series"""
    ## Set up the grid size
    #rows = 20  # 10x10 grid for 100 plots
    #cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    fig.suptitle(PATH_ANALYSIS, fontsize=16, y=0.995)  # Moved suptitle up
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)  # Adjust layout to leave space for suptitle

    for i in range(rows):
        for j in range(cols):
            index = coordid-1
            #print(f"Coordid {coordid}")

            if index < lambda_ref.shape[0]:
                ax = axes[i, j]
                
                if lambda_ref.iloc[index]['resname'] == "HSPT": #If histidine
                        if lambda_ref.iloc[index]['resid'] != prv_resid: #If this is histidine lambda1
                                
                            coordids = [coordid, coordid+1, coordid+2]
                            
                            res_prot_ts = get_histidine_protonation_timeseries(coordids, PATH_ANALYSIS, xvg_data)

                            prv_resid = lambda_ref.iloc[index]['resid'] 
                        else:
                            res_prot_ts = get_histidine_protonation_timeseries(coordids, PATH_ANALYSIS, xvg_data)

                else:
                    res_prot_ts = get_protonation_timeseries(coordid, PATH_ANALYSIS, xvg_data)
                
                if npz_output:
                    # Save output in npz file named after the residue
                    resname = lambda_ref.iloc[index]['resname']
                    resid = lambda_ref.iloc[index]['resid']
                    np.savez(f"{PATH_ANALYSIS}/{resname}_{resid}_protonation_timeseries.npz", res_prot_ts=res_prot_ts)

                ax.plot(res_prot_ts, label="MD1")
                ax.set_ylim(-0.1,1.1)

                # Set xticks and labels aaccording to the simulation time
                length = round(time/10000)*10000
                xticks = np.arange(0, length, length/3)
                xticks = np.round(xticks/10000)*10000
                xticks = np.concatenate((xticks, [length]))
                ax.set_xticks(xticks.astype(int))
                ax.set_xticklabels((xticks/1000).astype(int))
                
                ax.set_title(f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')

            else:
                continue
            
            coordid += 1

    for ax in axes.flat[coordid-1:]:
        ax.remove()

    
        
    return plt

def plot_protonation_convergence(PATH_ANALYSIS, time, xvg_data, coord2lambda_dict, lambda_ref, rows=20, cols=5):
    """Plot protonation avg and standard error time-series. Just two replicas supported - add get_protfrac_ts if more than 2"""
    # Set up the grid size
    #rows = 20  # 10x10 grid for 100 plots
    #cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    fig.suptitle(PATH_ANALYSIS, fontsize=16, y=0.995)  # Moved suptitle up
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)  # Adjust layout to leave space for suptitle

    for i in range(rows):
        for j in range(cols):
            index = coordid-1
            #print(f"Coordid {coordid}")

            if index < lambda_ref.shape[0]:
                ax = axes[i, j]
                
                if lambda_ref.iloc[index]['resname'] == "HSPT": #If histidine
                        if lambda_ref.iloc[index]['resid'] != prv_resid: #If this is histidine lambda1
                                
                            coordids = [coordid, coordid+1, coordid+2]
                            
                            res_prot_ts1 = get_histidine_protonation_timeseries(coordids, PATH_ANALYSIS[0], xvg_data)
                            res_prot_ts2 = get_histidine_protonation_timeseries(coordids, PATH_ANALYSIS[1], xvg_data)

                            prv_resid = lambda_ref.iloc[index]['resid'] 
                        else:
                            res_prot_ts1 = get_histidine_protonation_timeseries(coordids, PATH_ANALYSIS[0], xvg_data)
                            res_prot_ts2 = get_histidine_protonation_timeseries(coordids, PATH_ANALYSIS[1], xvg_data)
                else:
                    res_prot_ts1 = get_protonation_timeseries(coordid, PATH_ANALYSIS[0], xvg_data)
                    res_prot_ts2 = get_protonation_timeseries(coordid, PATH_ANALYSIS[1], xvg_data)
                
                min_length = min(len(res_prot_ts1), len(res_prot_ts2))
                total_protarray = np.vstack((res_prot_ts1[:min_length], res_prot_ts2[:min_length]))
                total_protonse = np.std(total_protarray, axis=0) / np.sqrt(len(total_protarray)) 
                ax.plot(np.mean(total_protarray, axis=0))
                ax.fill_between(np.arange(len(total_protarray[0,:])), np.mean(total_protarray, axis=0) - total_protonse, np.mean(total_protarray, axis=0) + total_protonse, alpha=0.5)

                # Set xticks and labels aaccording to the simulation time
                length = round(time/10000)*10000
                xticks = np.arange(0, length, length/3)
                xticks = np.round(xticks/10000)*10000
                xticks = np.concatenate((xticks, [length]))
                ax.set_xticks(xticks.astype(int))
                ax.set_xticklabels((xticks/1000).astype(int))

                ax.set_ylim(-0.1,1.1)
                ax.set_title(f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')

            else:
                continue
            
            coordid += 1

    for ax in axes.flat[coordid-1:]:
        ax.remove()

    return plt


def plot_protonation_fraction(PATH_ANALYSIS, xvg_data, lambda_ref, rows=20, cols=5, npz_output=False):
    """Plot protonation fractions avg and se"""
    # Set up the grid size
    #rows = 20  # 10x10 grid for 100 plots
    #cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout
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
            #print(f"Coordid {coordid}")

            if index < lambda_ref.shape[0]:
                ax = axes[i, j]
                
                # Generate some data; here, we're using sine waves with different frequencies

                if lambda_ref.iloc[index]['resname'] == "HSPT": #If histidine
                        if lambda_ref.iloc[index]['resid'] != prv_resid: #If this is histidine lambda1
                                
                            coordids = [coordid, coordid+1, coordid+2]
                            
                            proton_avg, proton_se = get_histidine_statistics(coordids, xvg_data)

                            prv_resid = lambda_ref.iloc[index]['resid'] 
                        else:
                            proton_avg, proton_se = get_histidine_statistics(coordids, xvg_data)

                else:
                    proton_avg, deproton_avg, proton_se, deproton_se  = get_statistics(coordid, xvg_data)


                ax.bar(['Protonated'], [proton_avg], yerr=proton_se)


                ax.set_ylim(0,1)
                ax.set_title(f'{lambda_ref.iloc[coordid-1]["resname"]}_{lambda_ref.iloc[coordid-1]["resid"]}')

                ax.set_xticks([])
                ax.set_yticks([])

                # Save output in npz file named after the residue
                if npz_output:
                    resname = lambda_ref.iloc[index]['resname']
                    resid = lambda_ref.iloc[index]['resid']
                    np.savez(f"{npz_dir}/{resname}_{resid}_protonation_fraction.npz", res_prot_avg=proton_avg, res_prot_se=proton_se)

            else:
                continue

            
            coordid += 1

    # Show the plot
    for ax in axes.flat[coordid-1:]:
        ax.remove()
        
    return plt

def single_residue_convergence(coordid, PATH_ANALYSIS, xvg_data, lambda_ref, title="Constant-pH MD"):
    """THIS WORKS ONLY FOR NON HISTIDINES! Plot convergence of single residue. Just two replicas supported - add res_fracX if more."""
    res_frac1 = get_protonation_timeseries(coordid, PATH_ANALYSIS[0], xvg_data)
    res_frac2 = get_protonation_timeseries(coordid, PATH_ANALYSIS[1], xvg_data)

    min_length = min(len(res_frac1), len(res_frac2))

    total_protarray = np.vstack((res_frac1[:min_length], res_frac2[:min_length]))

    total_protonse = np.std(total_protarray, axis=0) / np.sqrt(len(total_protarray)) 

    plt.plot(np.mean(total_protarray, axis=0), label=f'{lambda_ref.iloc[coordid-1]["resname"]}_{lambda_ref.iloc[coordid-1]["resid"]}')
    plt.fill_between(np.arange(len(total_protarray[0,:])), np.mean(total_protarray, axis=0) - total_protonse, np.mean(total_protarray, axis=0) + total_protonse, alpha=0.5)

    plt.ylim(0,1)
    plt.ylabel("Protonation fraction")
    plt.xlabel("Time (ps)")
    plt.legend()
    plt.title(title)
    return plt