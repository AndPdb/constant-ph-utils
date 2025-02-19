# If you want to run this outside of the constant-ph-utils directory, uncomment the following lines
import sys
# caution: path[0] is reserved for script path (or '' in REPL)
sys.path.insert(1, 'constant-ph-utils/')

import cProfile
import pstats
import argparse
from analyses import *
from plot import *
import pandas as pd

# Import lamdareference.dat

def main():
    lambda_ref = pd.read_csv('test/lambdareference.dat', sep=r'\s+', engine='python')

    coord2lambda_dict = parse_lambda_reference('test/lambdareference.dat')

    PATH_MD1 = "test/MD1/analysis"
    PATH_MD2 = "test/MD1_2/analysis"
    PATHS_MD = [PATH_MD1, PATH_MD2]
    MD1_PREFIX = "MD1"
    MD2_PREFIX = "MD1_2"

    # Create an instance of XVGData for the directories
    xvg_data = XVGData(PATHS_MD)

    # MD1

    ## Overview lambda distributions
    lambda_hist_md1 = plot_lambda_hist(PATH_MD1, xvg_data, coord2lambda_dict, lambda_ref)
    lambda_hist_md1.savefig(f"{MD1_PREFIX}_histograms.png")
    lambda_hist_md1.close()
    

    ## Protonation fraction time series
    proton_ts_md1 = plot_protonation_timeseries(PATH_MD1, xvg_data, coord2lambda_dict, lambda_ref)
    proton_ts_md1.savefig(f"{MD1_PREFIX}_timeseries.png")
    proton_ts_md1.close()

    # MD2

    # ## Overview lambda distributions
    lambda_hist_md2 = plot_lambda_hist(PATH_MD2, xvg_data, coord2lambda_dict, lambda_ref)
    lambda_hist_md2.savefig(f"{MD2_PREFIX}_histograms.png")
    lambda_hist_md2.close()

    # ## Protonation fraction time series
    proton_ts_md2 = plot_protonation_timeseries(PATH_MD2, xvg_data, coord2lambda_dict, lambda_ref)
    proton_ts_md2.savefig(f"{MD2_PREFIX}_timeseries.png")
    proton_ts_md2.close()

    # MD1 vs MD2

    # ## Protonation convergence
    proton_conv = plot_protonation_convergence(PATHS_MD, xvg_data, coord2lambda_dict, lambda_ref)
    proton_conv.savefig(f"{MD1_PREFIX}_{MD2_PREFIX}_convergence.png")
    proton_conv.close()

    # ## Overview protonation fractions
    proton_frac = plot_protonation_fraction(PATHS_MD, xvg_data, lambda_ref)
    proton_frac.savefig(f"{MD1_PREFIX}_{MD2_PREFIX}_protonfraction.png")
    proton_frac.close()

    # ## Sigle residue protonation fraction time series
    # ### Glu513
    res67_conv = single_residue_convergence(67, PATHS_MD, xvg_data, lambda_ref)
    res67_conv.savefig("Glu513.png")
    res67_conv.close()

    # ### Glu75
    res3_conv = single_residue_convergence(3, PATHS_MD, xvg_data, lambda_ref)
    res3_conv.savefig("Glu75.png")
    res3_conv.close()
    
    # ### Glu75
    res4_conv = single_residue_convergence(4, PATHS_MD, xvg_data, lambda_ref)
    res4_conv.savefig("Glu78.png")
    res4_conv.close()

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
