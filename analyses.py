from dataclasses import dataclass, field
import numpy as np
import pandas as pd
import os
import glob
from typing import Dict, Tuple, List
from concurrent.futures import ProcessPoolExecutor, as_completed


def load_file_for_pool(filepath: str, coord_id: int, num_rows: int) -> Tuple[int, np.ndarray]:
    data = np.empty((num_rows, 2), dtype=np.float64)
    index = 0
    with open(filepath, 'r') as f:
        for line in f:
            if line[0] in ('#', '@'):
                continue
            values = line.split()
            data[index, 0] = float(values[0])
            data[index, 1] = float(values[1])
            index += 1
            if index >= num_rows:
                break
    return (coord_id, data[:index, :])


@dataclass
class XVGData:
    """
    A class to handle the loading and processing of XVG data files.
    """

    analysis_dir: str = field(default_factory=str)
    coordids : List[int] = field(default_factory=list)
    num_rows: int = field(default_factory=2000000)
    num_threads: int = field(default_factory=2)
    data: Dict[int, np.ndarray] = field(default_factory=dict)

    def __post_init__(self):
        # safe place to call the loader after the instance exists
        self.load_all_data()

    def load_all_data(self):
        """
        Loads all XVG data from the specified directories into memory.
        """
        with ProcessPoolExecutor(max_workers=self.num_threads) as executor:
            futures_dict = {}

            #xvg_files = glob.glob(f"{self.analysis_dir}/*.xvg")
            # for analysis_dir in self.analysis_dirs:
            #for coord_id in range(1, len(xvg_files) + 1):
            for coord_id in self.coordids:
                coord_xvg_name = f"cphmd-coord-{coord_id}.xvg"
                coord_xvg_path = os.path.join(
                    self.analysis_dir, coord_xvg_name)
                print(coord_xvg_path)

                if os.path.exists(coord_xvg_path):
                    future = executor.submit(
                        load_file_for_pool, coord_xvg_path, coord_id, self.num_rows
                    )
                    futures_dict[future] = coord_id

            # Process results as they complete
            for future in as_completed(futures_dict):
                coord_id, arr = future.result()
                self.data[coord_id] = arr

    def __getitem__(self, coord_id: int | slice | tuple):

        if isinstance(coord_id, tuple):
            #coord_id = (row_index, col_id)
            row_index, col_id = coord_id
            column_data = self.data[col_id]
            return column_data[row_index] 
        
        elif isinstance(coord_id, slice):
            # Handle slicing
            start, stop, step = coord_id.indices(len(self.data))
            return [self.data[i] for i in range(start, stop, step)]
        
        elif isinstance(coord_id, int):
            # Handle single integer index
            return self.data[coord_id]
        
        else:
            raise TypeError(f"Invalid argument type: {type(coord_id)}")
    
    def __len__(self) -> int:
        """Returns the number of loaded coordinate files."""
        return len(self.data)

    def __contains__(self, coord_id: int) -> bool:
        """Check if a coordinate ID exists in the loaded data."""
        return coord_id in self.data


def calculate_fractions(array_xvg: np.ndarray) -> Tuple[float, float]:
    """
    Calculates protonation and deprotonation fractions from the given data.
    """
    # Identify protonated and deprotonated states based on a threshold
    prot = array_xvg[:, 1] < 0.2
    deprot = array_xvg[:, 1] > 0.8

    # Calculate total number of frames in protonated or deprotonated states
    total = prot.sum() + deprot.sum()

    # Handle the case where there are no protonated or deprotonated frames
    if total == 0:
        return 0.0, 0.0

    # Calculate fractions
    prot_frac = prot.sum() / total
    deprot_frac = deprot.sum() / total
    return prot_frac, deprot_frac


def get_statistics(coord_id: int, xvg_data_list: List[XVGData]) -> Tuple[float, float, float, float]:
    """
    Reads XVG files from replicas and computes statistics on protonation fractions.
    Returns averages and standard errors for protonation and deprotonation.
    """
    # Initialize lists to store fractions from different replicas
    prot_fractions = []
    deprot_fractions = []
    # Iterate over directories (replicas)
    for xvg_data in xvg_data_list:
        array_xvg = xvg_data[coord_id]
        if array_xvg.size > 0:
            prot_frac, deprot_frac = calculate_fractions(array_xvg)
            prot_fractions.append(prot_frac)
            deprot_fractions.append(deprot_frac)
    # Calculate average and standard error for protonation and deprotonation
    prot_avg = np.mean(prot_fractions)
    deprot_avg = np.mean(deprot_fractions)
    # ddof=1 use N-1 in the denominator
    prot_se = np.std(prot_fractions, ddof=1) / np.sqrt(len(prot_fractions))
    deprot_se = np.std(deprot_fractions, ddof=1) / \
        np.sqrt(len(deprot_fractions))
    return prot_avg, deprot_avg, prot_se, deprot_se


# def get_multichain_statistics(coord_id: int, xvg_data: XVGData) -> Tuple[float, float, float, float]:
#     """
#     Reads XVG files from replicas and computes statistics on protonation fractions.
#     Returns averages and standard errors for protonation and deprotonation.
#     """
#     # Initialize lists to store fractions from different replicas
#     prot_fractions = []
#     deprot_fractions = []

#     # Iterate over directories (replicas)
#     for analysis_dir in xvg_data.analysis_dirs:
#         for chain_id in chains:
#             array_xvg = xvg_data.get_coord_xvg(coord_id, analysis_dir)
#             if array_xvg.size > 0:
#                 prot_frac, deprot_frac = calculate_fractions(array_xvg)
#                 prot_fractions.append(prot_frac)
#                 deprot_fractions.append(deprot_frac)

#     # Calculate average and standard error for protonation and deprotonation
#     prot_avg = np.mean(prot_fractions)
#     deprot_avg = np.mean(deprot_fractions)
#     # ddof=1 use N-1 in the denominator
#     prot_se = np.std(prot_fractions, ddof=1) / np.sqrt(len(prot_fractions))
#     deprot_se = np.std(deprot_fractions, ddof=1) / \
#         np.sqrt(len(deprot_fractions))
#     return prot_avg, deprot_avg, prot_se, deprot_se


def get_protonation_timeseries(coord_id: int, xvg_data: XVGData) -> np.ndarray:
    """
    Computes the protonation fraction time series for the given coord_id.
    """
    # Read the XVG data
    array_xvg = xvg_data[coord_id]
    # Identify protonated and deprotonated states
    prot = array_xvg[:, 1] < 0.2
    deprot = array_xvg[:, 1] > 0.8
    # Calculate cumulative sums of protonated and deprotonated states
    prot_cumsum = np.cumsum(prot)
    deprot_cumsum = np.cumsum(deprot)
    total_cumsum = prot_cumsum + deprot_cumsum
    # Calculate the protonation fraction time series
    return prot_cumsum / total_cumsum


def calculate_histidine_fractions(array_xvgs: List[np.ndarray]) -> float:
    """
    Calculates the protonation fraction for histidine states.

    Args:
        array_xvgs: A list of three NumPy arrays, each representing a different histidine state (HSP, HSD, HSE).

    Returns:
        The protonation fraction of the histidine residue.
    """

    # Sum the number of neutral (nde) and protonated (npr) states
    nde = np.sum(array_xvgs[1][:, 1] > 0.8) + np.sum(array_xvgs[2][:, 1] > 0.8)
    npr = np.sum(array_xvgs[0][:, 1] > 0.8)

    # Calculate the protonation fraction
    prot_frac = npr / (nde + npr)

    return prot_frac


def get_histidine_statistics(coord_ids: List[int], xvg_data_list: List[XVGData]) -> Tuple[float, float]:
    """
    Reads XVG files for histidines and computes statistics on protonation fractions.
    """
    # Initialize a list to store protonation fractions
    prot_fractions = []
    # Iterate over directories (replicas)
    for xvg_data in xvg_data_list:
        # Read XVG data for each histidine coordinate
        histidine_data = [xvg_data[coord_id] for coord_id in coord_ids]
        # Calculate the protonation fraction for this replica
        prot_frac = calculate_histidine_fractions(histidine_data)
        prot_fractions.append(prot_frac)
    # Calculate average and standard error for protonation
    prot_avg = np.mean(prot_fractions)
    prot_se = np.std(prot_fractions, ddof=1) / np.sqrt(len(prot_fractions))
    return prot_avg, prot_se


def get_histidine_protonation_timeseries(coord_ids: List[int], xvg_data: XVGData) -> np.ndarray:
    """
    Computes the time series of protonation fractions for histidines.

    Args:
        coord_ids: A list of coordinate IDs for the histidine residues.
        xvg_data: The XVGData object containing the XVG data.

    Returns:
        A NumPy array representing the protonation fraction time series for the histidine residue.
    """

    # Read XVG data for each histidine coordinate
    xvg_his_dict = {cid: xvg_data[cid] for cid in coord_ids}

    # Identify protonated and deprotonated states for each time step
    # If the first coordinate (HSP) is > 0.8, it's protonated
    protonated_states = xvg_his_dict[coord_ids[0]][:, 1] > 0.8
    # If either HSD or HSE is > 0.8, it's deprotonated
    deprotonated_states = np.any(
        [xvg_his_dict[cid][:, 1] > 0.8 for cid in coord_ids[1:]], axis=0)

    # Calculate cumulative sums of protonated and deprotonated states over time
    prot_cumsum = np.cumsum(protonated_states)
    deprot_cumsum = np.cumsum(deprotonated_states)

    # Calculate the protonation fraction time series, handling potential division by zero
    protonation_timeseries = prot_cumsum / (prot_cumsum + deprot_cumsum)

    return protonation_timeseries


def parse_lambda_reference(lambda_reference_path: str) -> Dict[int, int]:
    """
    Parses the lambda reference file and returns a dictionary mapping coord_id to lambda_id.
    """
    # Initialize a dictionary to store the mapping
    coord_to_lambda = {}
    # Read the lambda reference file
    with open(lambda_reference_path, "r") as file:
        coord_id = 0
        lambda_id = 0
        for line in file:
            # Skip header lines
            if line.startswith("resname"):
                continue
            # Increment coord_id and lambda_id. Count HIS states only once.
            if line[5] != "1":
                coord_id += 1
            else:
                lambda_id += 1
                coord_id += 1
            # Add the mapping to the dictionary
            coord_to_lambda[coord_id] = lambda_id
    return coord_to_lambda


def resid2coordid(resid: int, lambda_ref: pd.DataFrame) -> int:
    """
    Converts a residue ID to a coordinate ID using the lambda reference DataFrame.
    """

    # Extract the row where 'resid' is 75
    row = lambda_ref[lambda_ref['resid'] == resid]

    # Get the 'coordinateFile' value from the row
    coordinate_file = row['coordinateFile'].values[0]

    # Split the filename to extract the desired part
    file_parts = coordinate_file.split('-')
    coordinate_id = file_parts[-1].split('.xvg')[0]

    return int(coordinate_id)
