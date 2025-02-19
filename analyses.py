import numpy as np
import os, glob
from typing import Dict, Tuple, List
from concurrent.futures import ThreadPoolExecutor

class XVGData:
    """
    A class to handle the loading and processing of XVG data files from multiple directories.
    """

    def __init__(self, analysis_dirs: List[str], num_rows: int = 2000000):
        """
        Initializes the XVGData object.

        Args:
            analysis_dirs (List[str]): List of directories containing the XVG files.
            num_rows (int): Maximum number of rows to read from each XVG file.
        """
        self.analysis_dirs = analysis_dirs
        self.num_rows = num_rows
        self.data = {dir: {} for dir in analysis_dirs}
        self.load_all_data()

    def load_all_data(self):
        """
        Loads all XVG data from the specified directories into memory.
        """
        with ThreadPoolExecutor() as executor:
            futures = []
            for analysis_dir in self.analysis_dirs:
                for coord_id in range(1, len(glob.glob(f"{analysis_dir}/*.xvg")) + 1):
                    coord_xvg_name = f"cphmd-coord-{coord_id}.xvg"
                    coord_xvg_path = os.path.join(analysis_dir, coord_xvg_name)
                    print(coord_xvg_path)
                    if os.path.exists(coord_xvg_path):
                        futures.append(executor.submit(self._load_file, coord_xvg_path, analysis_dir, coord_id))
            for future in futures:
                future.result()

    def _load_file(self, filepath: str, analysis_dir: str, coord_id: int):
        """
        Helper function to load a single XVG file.
        """
        self.data[analysis_dir][coord_id] = self.fast_read_numpy_array(filepath)

    def fast_read_numpy_array(self, filename: str) -> np.ndarray:
        """
        Reads an XVG file into a NumPy array, skipping comment lines.

        Args:
            filename (str): Path to the XVG file.

        Returns:
            np.ndarray: Array containing the data from the XVG file.
        """
        data = np.empty((self.num_rows, 2), dtype=np.float64)  # Preallocate array
        index = 0  # Keep track of valid data row index

        with open(filename, 'r') as f:
            for line in f:
                if line[0] in ('#', '@'):  # Skip comments
                    continue
                values = line.split()
                data[index, 0] = float(values[0])
                data[index, 1] = float(values[1])
                index += 1
                if index >= self.num_rows:  # Prevent overflow
                    break

        return data[:index, :]
    
    def get_coord_xvg(self, coord_id: int, analysis_dir: str) -> np.ndarray:
        return self.data[analysis_dir].get(coord_id, np.array([]))
    
    
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


def get_statistics(coord_id: int, xvg_data: XVGData) -> Tuple[float, float, float, float]:
    """
    Reads XVG files from replicas and computes statistics on protonation fractions.
    Returns averages and standard errors for protonation and deprotonation.
    """
    # Initialize lists to store fractions from different replicas
    prot_fractions = []
    deprot_fractions = []
    # Iterate over directories (replicas)
    for analysis_dir in xvg_data.analysis_dirs:
        array_xvg = xvg_data.get_coord_xvg(coord_id, analysis_dir)
        if array_xvg.size > 0:
            prot_frac, deprot_frac = calculate_fractions(array_xvg)
            prot_fractions.append(prot_frac)
            deprot_fractions.append(deprot_frac)
    # Calculate average and standard error for protonation and deprotonation
    prot_avg = np.mean(prot_fractions)
    deprot_avg = np.mean(deprot_fractions)
    prot_se = np.std(prot_fractions, ddof=1) / np.sqrt(len(prot_fractions)) #ddof=1 use N-1 in the denominator 
    deprot_se = np.std(deprot_fractions, ddof=1) / np.sqrt(len(deprot_fractions))
    return prot_avg, deprot_avg, prot_se, deprot_se


def get_protonation_timeseries(coord_id: int, analysis_dir: str, xvg_data: XVGData) -> np.ndarray:
    """
    Computes the protonation fraction time series for the given coord_id.
    """
    # Read the XVG data
    array_xvg = xvg_data.get_coord_xvg(coord_id, analysis_dir)
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


def get_histidine_statistics(coord_ids: List[int], xvg_data: XVGData) -> Tuple[float, float]:
    """
    Reads XVG files for histidines and computes statistics on protonation fractions.
    """
    # Initialize a list to store protonation fractions
    prot_fractions = []
    # Iterate over directories (replicas)
    for analysis_dir in xvg_data.analysis_dirs:
        # Read XVG data for each histidine coordinate
        histidine_data = [xvg_data.get_coord_xvg(coord_id, analysis_dir) for coord_id in coord_ids]
        # Calculate the protonation fraction for this replica
        prot_frac = calculate_histidine_fractions(histidine_data)
        prot_fractions.append(prot_frac)
    # Calculate average and standard error for protonation
    prot_avg = np.mean(prot_fractions)
    prot_se = np.std(prot_fractions, ddof=1) / np.sqrt(len(prot_fractions))
    return prot_avg, prot_se


def get_histidine_protonation_timeseries(coord_ids: List[int], analysis_dir: str, xvg_data: XVGData) -> np.ndarray:
    """
    Computes the time series of protonation fractions for histidines.

    Args:
        coord_ids: A list of coordinate IDs for the histidine residues.
        analysis_dir: The directory containing the XVG files.

    Returns:
        A NumPy array representing the protonation fraction time series for the histidine residue.
    """

    # Read XVG data for each histidine coordinate
    xvg_his_dict = {cid: xvg_data.get_coord_xvg(cid, analysis_dir) for cid in coord_ids}

    # Identify protonated and deprotonated states for each time step
    protonated_states = xvg_his_dict[coord_ids[0]][:, 1] > 0.8  # If the first coordinate (HSP) is > 0.8, it's protonated
    deprotonated_states = np.any([xvg_his_dict[cid][:, 1] > 0.8 for cid in coord_ids[1:]], axis=0)  # If either HSD or HSE is > 0.8, it's deprotonated

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

def resid2coordid(resid: int, lambda_ref: Dict[int, int]) -> int:
    """
    Converts a residue ID to a coordinate ID using the lambda reference dictionary.
    """
    
    # Extract the row where 'resid' is 75
    row = lambda_ref[lambda_ref['resid'] == resid]

    # Get the 'coordinateFile' value from the row
    coordinate_file = row['coordinateFile'].values[0]

    # Split the filename to extract the desired part
    file_parts = coordinate_file.split('-')
    coordinate_id = file_parts[-1].split('.xvg')[0]

    return int(coordinate_id)
