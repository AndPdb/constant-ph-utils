# If you want to run this outside of the constant-ph-utils directory, uncomment the following lines
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, 'constant-ph-utils/')

import cProfile
import pstats
import argparse
from analyses import *
from plot import *
import matplotlib.pyplot as plt

# Import lamdareference.dat
import pandas as pd

def main():
    lambda_ref = pd.read_csv('test/lambdareference.dat', sep='\s+')

    coord2lambda_dict = parse_lambda_reference('test/lambdareference.dat')

    PATH_MD1 = "test/MD1/analysis"
    PATH_MD2 = "test/MD1_2/analysis"
    PATHS_MD = [PATH_MD1, PATH_MD2]
    MD1_PREFIX = "MD1"
    MD2_PREFIX = "MD1_2"

    # # MD1

    # ## Overview lambda distributions

    # Set up the grid size
    rows = 20  # 10x10 grid for 100 plots
    cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    fig.suptitle(PATH_MD1, fontsize=16, y=0.995)  # Moved suptitle up
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)  # Adjust layout to leave space for suptitle

    for i in range(rows):
        for j in range(cols):
            index = coordid-1
            #print(f"Coordid {coordid}")

            if index < lambda_ref.shape[0]:
                ax = axes[i, j]
                data = read_coord_xvg(str(coordid), PATH_MD1)
                ax.hist(data[:,1], bins=1000);

                ax.set_xlim(-0.125,1.125)
                ax.set_title(f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')

            else:
                continue
            
            coordid += 1

    # Show the plot
    #plt.show()
    plt.savefig(f"{MD1_PREFIX}_histograms.png")
    plt.close()

    # ## Protonation fraction time series


    # Set up the grid size
    rows = 20  # 10x10 grid for 100 plots
    cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    fig.suptitle(PATH_MD1, fontsize=16, y=0.995)  # Moved suptitle up
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
                            
                            res_prot_ts = get_histidine_protonation_timeseries(coordids, PATH_MD1)

                            prv_resid = lambda_ref.iloc[index]['resid'] 
                        else:
                            res_prot_ts = get_histidine_protonation_timeseries(coordids, PATH_MD1)

                else:
                    res_prot_ts = get_protonation_timeseries(coordid, PATH_MD1)

                ax.plot(res_prot_ts, label="MD1")
                ax.set_ylim(-0.1,1.1)
                ax.set_title(f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')

            else:
                continue
            
            coordid += 1

    # Show the plot
    #plt.show()
    plt.savefig(f"{MD1_PREFIX}_timeseries.png")
    plt.close()

    # # MD1_3

    # ## Overview lambda distributions

    # Set up the grid size
    rows = 20  # 10x10 grid for 100 plots
    cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    fig.suptitle(PATH_MD2, fontsize=16, y=0.995)  # Moved suptitle up
    fig.tight_layout(rect=[0, 0, 1, 0.99], pad=2.0)  # Adjust layout to leave space for suptitle

    for i in range(rows):
        for j in range(cols):
            index = coordid-1
            #print(f"Coordid {coordid}")

            if index < lambda_ref.shape[0]:
                ax = axes[i, j]
                data = read_coord_xvg(str(coordid), PATH_MD2)
                ax.hist(data[:,1], bins=1000);
                
                ax.set_xlim(-0.125,1.125)
                ax.set_title(f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')

            else:
                continue
            
            coordid += 1

    # Show the plot
    #plt.show()
    plt.savefig(f"{MD2_PREFIX}_histograms.png")
    plt.close()

    # ## Protonation fraction time series

    # Set up the grid size
    rows = 20  # 10x10 grid for 100 plots
    cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    fig.suptitle(PATH_MD2, fontsize=16, y=0.995)  # Moved suptitle up
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
                            
                            res_prot_ts = get_histidine_protonation_timeseries(coordids, PATH_MD2)

                            prv_resid = lambda_ref.iloc[index]['resid'] 
                        else:
                            res_prot_ts = get_histidine_protonation_timeseries(coordids, PATH_MD2)

                else:
                    res_prot_ts = get_protonation_timeseries(coordid, PATH_MD2)

                ax.plot(res_prot_ts, label="MD1_3")
                ax.set_ylim(-0.1,1.1)
                ax.set_title(f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')

            else:
                continue
            
            coordid += 1

    # Show the plot
    #plt.show()
    plt.savefig(f"{MD2_PREFIX}_timeseries.png")
    plt.close()

    # MD1 vs MD1_3


    # Set up the grid size
    rows = 20  # 10x10 grid for 100 plots
    cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout

    coordid = 1
    prv_resid = 0

    # Generate and plot data for each subplot
    fig.suptitle(PATHS_MD, fontsize=16, y=0.995)  # Moved suptitle up
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
                            
                            res_prot_ts1 = get_histidine_protonation_timeseries(coordids, PATHS_MD[0])
                            res_prot_ts2 = get_histidine_protonation_timeseries(coordids, PATHS_MD[1])

                            prv_resid = lambda_ref.iloc[index]['resid'] 
                        else:
                            res_prot_ts1 = get_histidine_protonation_timeseries(coordids, PATHS_MD[0])
                            res_prot_ts2 = get_histidine_protonation_timeseries(coordids, PATHS_MD[1])
                else:
                    res_prot_ts1 = get_protonation_timeseries(coordid, PATHS_MD[0])
                    res_prot_ts2 = get_protonation_timeseries(coordid, PATHS_MD[1])
                
                min_length = min(len(res_prot_ts1), len(res_prot_ts2))
                total_protarray = np.vstack((res_prot_ts1[:min_length], res_prot_ts2[:min_length]))
                total_protonse = np.std(total_protarray, axis=0) / np.sqrt(len(total_protarray)) 
                ax.plot(np.mean(total_protarray, axis=0))
                ax.fill_between(np.arange(len(total_protarray[0,:])), np.mean(total_protarray, axis=0) - total_protonse, np.mean(total_protarray, axis=0) + total_protonse, alpha=0.5)
                ax.set_ylim(-0.1,1.1)
                ax.set_title(f'coord_{coordid}-lambda_{coord2lambda_dict[coordid]}')

            else:
                continue
            
            coordid += 1

    # Show the plot
    #plt.show()
    plt.savefig(f"{MD1_PREFIX}_{MD2_PREFIX}_convergence.png")
    plt.close()


    # ## Overview protonation fractions

    # Set up the grid size
    rows = 20  # 10x10 grid for 100 plots
    cols = 5

    # Create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(15, 40))  # Adjusted figure size for 5x20 layout
    fig.tight_layout(pad=2.0)  # Optional: adjust padding between plots

    coordid = 1
    prv_resid = 0



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
                            
                            proton_avg, proton_se = get_histidine_statistics(coordids, PATHS_MD)

                            prv_resid = lambda_ref.iloc[index]['resid'] 
                        else:
                            proton_avg, proton_se = get_histidine_statistics(coordids, PATHS_MD)

                else:
                    proton_avg, deproton_avg, proton_se, deproton_se  = get_statistics(str(coordid), PATHS_MD)


                ax.bar(['Protonated'], [proton_avg], yerr=proton_se)


                ax.set_ylim(0,1)
                ax.set_title(f'{lambda_ref.iloc[coordid-1]["resname"]}_{lambda_ref.iloc[coordid-1]["resid"]}')

                ax.set_xticks([])
                ax.set_yticks([])
            else:
                continue

            
            coordid += 1

    # Show the plot
    #plt.show()
    plt.savefig(f"{MD1_PREFIX}_{MD2_PREFIX}_protonfraction.png")
    plt.close()

    # ## Sigle residue protonation fraction time series

    # ### Glu513
    glu513_frac1 = get_protonation_timeseries("67", PATHS_MD[0])
    glu513_frac2 = get_protonation_timeseries("67", PATHS_MD[1])

    min_length = min(len(glu513_frac1), len(glu513_frac2))

    total_protarray = np.vstack((glu513_frac1[:min_length], glu513_frac2[:min_length]))

    total_protonse = np.std(total_protarray, axis=0) / np.sqrt(len(total_protarray)) 

    plt.plot(np.mean(total_protarray, axis=0), label="Glu513")
    plt.fill_between(np.arange(len(total_protarray[0,:])), np.mean(total_protarray, axis=0) - total_protonse, np.mean(total_protarray, axis=0) + total_protonse, alpha=0.5)

    plt.ylim(0,1)
    plt.ylabel("Protonation fraction")
    plt.xlabel("Time (ps)")
    plt.legend()
    plt.title("GTR1 IF cphmd at pH 7.5")
    #plt.show()
    plt.savefig("Glu513.png")
    plt.close()


    # ### Glu75
    glu75_frac1 = get_protonation_timeseries("3", PATHS_MD[0])
    glu75_frac2 = get_protonation_timeseries("3", PATHS_MD[1])

    min_length = min(len(glu75_frac1), len(glu75_frac2))

    total_protarray = np.vstack((glu75_frac1[:min_length], glu75_frac2[:min_length]))

    total_protonse = np.std(total_protarray, axis=0) / np.sqrt(len(total_protarray)) 

    plt.plot(np.mean(total_protarray, axis=0), label="Glu75")
    plt.fill_between(np.arange(len(total_protarray[0,:])), np.mean(total_protarray, axis=0) - total_protonse, np.mean(total_protarray, axis=0) + total_protonse, alpha=0.5)

    plt.ylim(0,1)
    plt.ylabel("Protonation fraction")
    plt.xlabel("Time (ps)")
    plt.legend()
    plt.title("GTR1 IF cphmd at pH 7.5")
    #plt.show()
    plt.savefig("Glu75.png")
    plt.close()

    # ### Glu78
    glu78_frac1 = get_protonation_timeseries("4", PATHS_MD[0])
    glu78_frac2 = get_protonation_timeseries("4", PATHS_MD[1])

    min_length = min(len(glu78_frac1), len(glu78_frac2))

    total_protarray = np.vstack((glu78_frac1[:min_length], glu78_frac2[:min_length]))

    total_protonse = np.std(total_protarray, axis=0) / np.sqrt(len(total_protarray)) 

    plt.plot(np.mean(total_protarray, axis=0), label="Glu78")
    plt.fill_between(np.arange(len(total_protarray[0,:])), np.mean(total_protarray, axis=0) - total_protonse, np.mean(total_protarray, axis=0) + total_protonse, alpha=0.5)

    plt.ylim(0,1)
    plt.ylabel("Protonation fraction")
    plt.xlabel("Time (ps)")
    plt.legend()
    plt.title("GTR1 IF cphmd at pH 7.5")
    #plt.show()
    plt.savefig("Glu78.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run analysis script with optional profiling.")
    parser.add_argument('--profile', action='store_true', help="Enable profiling with cProfile")
    args = parser.parse_args()

    if args.profile:
        profiler = cProfile.Profile()
        profiler.enable()
        main()
        profiler.disable()
        profiler.dump_stats('profile_results.prof')

        # Optionally, print the profiling results to a text file
        with open('profile_results.txt', 'w') as f:
            stats = pstats.Stats(profiler, stream=f)
            stats.sort_stats('cumulative')
            stats.print_stats()
    else:
        main()