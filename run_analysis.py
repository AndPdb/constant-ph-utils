# PYTHON_ARGCOMPLETE_OK

import sys, os
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, 'constant-ph-utils/')

import pandas as pd
from plot import *
from analyses import *
import argparse
import argcomplete
import pstats
import cProfile


# RUN_TYPE = "Publication"  # Debug or Publication
# PLOT_TYPE = "Debug"  # Debug or Publication
# SINGLE_LETTER = True  # Use single-letter amino acid codes in labels
# PLOT_ROWS = 21  # Number of rows in the plot grid
# PLOT_COLS = 5  # Number of columns in the plot grid

# NPZ_OUTPUT = False  # True or False
# THREADS = 8  # Number of threads for parallel processing
# XVG_ROWS = 2000000  # Number of rows to read from the XVG files

# # List of chains to analyze. Right now works for two chains only.
# CHAINS = None

# # Path variables
# LAMBDAREF_PATH = "test"
# PATH_MD1 = "test/MD1/analysis"
# PATH_MD2 = "test/MD1_2/analysis"
# PATHS_MD = [PATH_MD1, PATH_MD2]
# MD1_PREFIX = "MD1"
# MD2_PREFIX = "MD1_2"
# res_ids = [75, 78, 513] # Residue number for single convergence plot

# CONVERG_PREFIX = "-".join([x.split("/")[-2] for x in PATHS_MD])


###### Main function #######


def main(args):
    # Unpack arguments into local variables
    RUN_TYPE = args.run_type
    PLOT_TYPE = args.plot_type
    OUTPUT_DIR_PLOT = args.dir_plot
    RES_IDS = args.res_ids
    SINGLE_LETTER = args.single_letter
    PLOT_ROWS = args.plot_rows
    PLOT_COLS = args.plot_cols
    NPZ_OUTPUT = args.npz_output
    THREADS = args.threads
    XVG_ROWS = args.xvg_rows
    CHAINS = args.chains
    LAMBDAREF_PATH = args.lambdaref_path
    PATHS_MD = args.paths_md

 
    # Derive convergence prefix from MD paths
    CONVERG_PREFIX = "-".join([x.rstrip("/").split("/")[-2] for x in PATHS_MD])
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
        title = path.rstrip("/").split("/")[-2]
        if CHAINS is not None:
            for chain in lambda_ref_chains.groups.keys():
                lambda_hist = plot_lambda_hist(
                    xvg_data_list[i], coord2lambda_dict, lambda_ref, rows=PLOT_ROWS, cols=PLOT_COLS, single_letter=SINGLE_LETTER)
                lambda_hist.savefig(os.path.join(OUTPUT_DIR_PLOT,f"{title}_{chain}_histograms.png"))
                lambda_hist.close()
                i += 1
        else:
            lambda_hist = plot_lambda_hist(
                xvg_data_list[i], coord2lambda_dict, lambda_ref, rows=PLOT_ROWS, cols=PLOT_COLS,  single_letter=SINGLE_LETTER)
            lambda_hist.savefig(os.path.join(OUTPUT_DIR_PLOT,f"{title}_histograms.png"))
            lambda_hist.close()
            i += 1

    # Protonation fraction time series
    i = 0
    for path in PATHS_MD:
        title = path.rstrip("/").split("/")[-2]
        if CHAINS is not None:
            for chain in CHAINS:
                proton_ts = plot_protonation_timeseries(
                    time_MD1, xvg_data_list[i], coord2lambda_dict, lambda_ref, rows=PLOT_ROWS, single_letter=SINGLE_LETTER)
                proton_ts.savefig(os.path.join(OUTPUT_DIR_PLOT,f"{title}_{chain}_timeseries.png"))
                proton_ts.close()
                i += 1
        else:
            proton_ts = plot_protonation_timeseries(
                time_MD1, xvg_data_list[i], coord2lambda_dict, lambda_ref, rows=PLOT_ROWS, single_letter=SINGLE_LETTER)
            proton_ts.savefig(os.path.join(OUTPUT_DIR_PLOT,f"{title}_timeseries.png"))
            proton_ts.close()
            i += 1

    if RUN_TYPE == "Publication":     # In case we have replicas

        # ## Protonation convergence
        proton_conv = plot_protonation_convergence(
            PATHS_MD, min_time, xvg_data_list, coord2lambda_dict, lambda_ref, chain_mapping=mapping, rows=PLOT_ROWS, quality=PLOT_TYPE,  single_letter=SINGLE_LETTER)
        proton_conv.savefig(os.path.join(OUTPUT_DIR_PLOT,f"{CONVERG_PREFIX}_convergence.png"), dpi=300)
        proton_conv.close()

        # ## Overview protonation fractions
        fig = plot_protonation_fraction(
            xvg_data_list, lambda_ref, chain_mapping=mapping,
            npz_output=NPZ_OUTPUT,  single_letter=SINGLE_LETTER)
        fig.savefig(os.path.join(OUTPUT_DIR_PLOT,f"{CONVERG_PREFIX}_protonfraction.png"), bbox_inches='tight', dpi=300)
        plt.close(fig)

        # ## Sigle residue protonation fraction time series
        for res_id in RES_IDS:
            res_coord = resid2coordid(res_id, lambda_ref)
            res_conv = single_residue_convergence(
                res_coord, xvg_data_list, lambda_ref, chain_mapping=mapping, single_letter=SINGLE_LETTER)
            res_conv.savefig(os.path.join(OUTPUT_DIR_PLOT, f"Res_{res_id}.png"),dpi=300)
            res_conv.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run constant-pH MD analysis and generate plots.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
 
    # --- Paths (required) ---
    parser.add_argument('--lambdaref-path', type=str, required=True,
                        help="Directory containing lambdareference.dat")
    parser.add_argument('--paths-md', type=str, nargs='+', required=True,
                        help="One or more paths to MD analysis directories")
    parser.add_argument('--dir-plot', type=str, required=True,
                        help="Output folder for plots")
    # --- Single-residue convergence ---
    parser.add_argument('--res-ids', nargs='+', type=int, default=None,
                        help="Residue IDs for single-residue convergence plots")

    # --- Run / plot mode ---
    parser.add_argument('--run-type', choices=['Debug', 'Publication'],
                        default='Publication',
                        help="Debug runs single-replica plots only; "
                             "Publication adds convergence and fraction plots")
    parser.add_argument('--plot-type', choices=['Debug', 'Publication'],
                        default='Debug',
                        help="Debug shows coord IDs in titles; "
                             "Publication shows residue names")
 
    # --- Grid dimensions ---
    parser.add_argument('--plot-rows', type=int, default=21,
                        help="Number of rows in the overview plot grids")
    parser.add_argument('--plot-cols', type=int, default=5,
                        help="Number of columns in the overview plot grids")
 
    # --- Labels ---
    parser.add_argument('--no-single-letter', dest='single_letter',
                        action='store_false',
                        help="Use three-letter amino acid codes instead of one-letter")
    parser.set_defaults(single_letter=True)
 
    # --- Output options ---
    parser.add_argument('--npz-output', action='store_true', default=False,
                        help="Save protonation data as .npz files")
    parser.add_argument('--dpi', type=int, default=300,
                        help="DPI for saved figures")
 
    # --- Performance ---
    parser.add_argument('--threads', type=int, default=8,
                        help="Number of threads for parallel XVG loading")
    parser.add_argument('--xvg-rows', type=int, default=2000000,
                        help="Maximum number of rows to read per XVG file")
 
    # --- Multi-chain ---
    parser.add_argument('--chains', nargs='+', default=None,
                        help="Chain identifiers for homomeric systems (e.g. A B)")
 

 
    # --- Profiling ---
    parser.add_argument('--profile', action='store_true',
                        help="Enable profiling with cProfile")

    argcomplete.autocomplete(parser) 
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
        main(args)
