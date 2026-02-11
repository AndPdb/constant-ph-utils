# If you want to run this outside of the constant-ph-utils directory, uncomment the following lines
import sys, os
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, 'constant-ph-utils/')

import pandas as pd
from plot import *
from analyses import *
import argparse
import pstats
import cProfile


RUN_TYPE = "Publication"  # Debug or Publication
PLOT_TYPE = "Debug"  # Debug or Publication
NPZ_OUTPUT = False  # True or False
THREADS = 8  # Number of threads for parallel processing
XVG_ROWS = 2000000  # Number of rows to read from the XVG files
PLOT_ROWS = 21  # Number of rows in the plot grid
PLOT_COLS = 5  # Number of columns in the plot grid
# List of chains to analyze. Right now works for two chains only.
CHAINS = None

# Path variables
LAMBDAREF_PATH = "test"
PATH_MD1 = "test/MD1/analysis"
PATH_MD2 = "test/MD1_2/analysis"
PATHS_MD = [PATH_MD1, PATH_MD2]
MD1_PREFIX = "MD1"
MD2_PREFIX = "MD1_2"

###### Main function #######


def main():
    # Load lambda reference data
    lambda_ref = pd.read_csv(os.path.join(
        LAMBDAREF_PATH, 'lambdareference.dat'), sep=r'\s+', engine='python')

    if lambda_ref['chain'][0].endswith("xvg"):
        lambda_ref = lambda_ref.rename(
            columns={"chain": "coordinateFile", "coordinateFile": "chain"})

    # Create a mapping from coordinate IDs to lambda values
    coord2lambda_dict = parse_lambda_reference(
        os.path.join(LAMBDAREF_PATH, 'lambdareference.dat'))

    # Get the list of coordinate IDs
    coordids = list(coord2lambda_dict.keys())

    # Create an instance of XVGData for the directories
    xvg_data_list = []

    # In case of multiple homomeric chains
    if CHAINS is not None:
        lambda_ref_chains = lambda_ref.groupby('chain')
        # This should be adapted case by case, depending on the chain names in the lambda reference and the desired mapping
        mapping = dict(
            zip(lambda_ref_chains.groups[CHAINS[0]]+1, lambda_ref_chains.groups[CHAINS[1]]+1))

        coordids_chain = {}
        for chain in CHAINS:
            coordids_chain[str(chain)] = [
                index+1 for index, row in lambda_ref_chains.get_group(str(chain)).iterrows()]

    # In case of single chain
    else:
        mapping = {}

    # Load XVGdata for each path and chain
    for path in PATHS_MD:
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"Directory {path} does not exist. Please check the path and try again.")

        if CHAINS is not None:
            for chain in CHAINS:
                xvg_data_list.append(XVGData(path, coordids=coordids_chain[str(
                    chain)], num_rows=XVG_ROWS, num_threads=THREADS))
        else:
            xvg_data_list.append(
                XVGData(path, coordids=coordids, num_rows=XVG_ROWS, num_threads=THREADS))

    # Get last time of MD1 for plotting time series
    time_MD1 = xvg_data_list[0][1][-1, 0]

    # In case we have replicas
    if RUN_TYPE == "Publication":
        time_MDs = []

        for data in xvg_data_list:
            try:
                time_MDs.append(data[1][-1, 0])
            except KeyError:
                continue
        # Find common minimal time across all simulations for convergence analysis
        min_time = min(time_MDs)

    # Overview lambda distributions
    i = 0
    for path in PATHS_MD:
        title = path.split("/")[-2]
        if CHAINS is not None:
            for chain in lambda_ref_chains.groups.keys():
                lambda_hist = plot_lambda_hist(
                    xvg_data_list[i], coord2lambda_dict, lambda_ref, rows=PLOT_ROWS, cols=PLOT_COLS)
                lambda_hist.savefig(f"{title}_{chain}_histograms.png")
                lambda_hist.close()
                i += 1
        else:
            lambda_hist = plot_lambda_hist(
                xvg_data_list[i], coord2lambda_dict, lambda_ref, rows=PLOT_ROWS, cols=PLOT_COLS)
            lambda_hist.savefig(f"{title}_histograms.png")
            lambda_hist.close()
            i += 1

    # Protonation fraction time series
    i = 0
    for path in PATHS_MD:
        title = path.split("/")[-2]
        if CHAINS is not None:
            for chain in CHAINS:
                proton_ts = plot_protonation_timeseries(
                    time_MD1, xvg_data_list[i], coord2lambda_dict, lambda_ref, rows=PLOT_ROWS)
                proton_ts.savefig(f"{title}_{chain}_timeseries.png")
                proton_ts.close()
                i += 1
        else:
            proton_ts = plot_protonation_timeseries(
                time_MD1, xvg_data_list[i], coord2lambda_dict, lambda_ref, rows=PLOT_ROWS)
            proton_ts.savefig(f"{title}_timeseries.png")
            proton_ts.close()
            i += 1

    if RUN_TYPE == "Publication":     # In case we have replicas

        # ## Protonation convergence
        proton_conv = plot_protonation_convergence(
            PATHS_MD, min_time, xvg_data_list, coord2lambda_dict, lambda_ref, chain_mapping=mapping, rows=PLOT_ROWS, quality=PLOT_TYPE)
        proton_conv.savefig(f"{MD1_PREFIX}_{MD2_PREFIX}_convergence.png")
        proton_conv.close()

        # ## Overview protonation fractions
        proton_frac = plot_protonation_fraction(
            xvg_data_list, lambda_ref, chain_mapping=mapping, rows=PLOT_ROWS, npz_output=NPZ_OUTPUT)
        proton_frac.savefig(f"{MD1_PREFIX}_{MD2_PREFIX}_protonfraction.png")
        proton_frac.close()

        # ## Sigle residue protonation fraction time series
        # ### Glu513
        glu513_id = resid2coordid(513, lambda_ref)
        res3_conv = single_residue_convergence(
            glu513_id, xvg_data_list, lambda_ref, chain_mapping=mapping, title="Convergence of residue Glu513")
        res3_conv.savefig("Glu513.png")
        res3_conv.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run analysis script with optional profiling.")
    parser.add_argument('--profile', action='store_true',
                        help="Enable profiling with cProfile")
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
