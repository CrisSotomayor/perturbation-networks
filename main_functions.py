"""Python scripts with the main functions used across the jupyter notebooks."""
#!/bin/user/env python3


import os
import numpy as np
import scipy as sp
import pandas as pd
import biographs as bg
import networkx as nx
import Bio.PDB.Polypeptide as pp

from collections import Counter


# Some global variables
AA = list(pp.aa1)
proteins = ['1be9', '1d5r', '1nd4', '3dqw', '4bz3']
protein_names = ['PSD95', 'PTEN', 'APH(3\')II', 'Src CD', 'VIM-2']
measures = ['nodes', 'edges', 'weight', 'distance']


# Get functional data
# Import processed functional data as DataFrames, all files have ordered AA list
# as index, positions as columns. Save data in functional_data.
functional_data = dict()
for protein in proteins:
    csv_file = os.path.join('data', f'functional_{protein}.csv')
    functional_data[protein] = pd.read_csv(csv_file, index_col=0, header=0)


def ReadNetworkCSV(protein, threshold, measure, data_path='data/structure'):
    """Return DataFrame from corresponding CSV.

    If protein has multiple identical chains, return average value
    for each position amongst all chains.


    Arguments
    --------
    protein: str, protein PDB name.
    threshold: float, distance threshold in Angstroms.
    measure: str, one of the structural measures.
    data_path (keyword): str, path to directory with pdbs.
    """

    file = os.path.join(data_path,
                        f"{protein}/{protein}_{threshold}_{measure}.csv")
    network_df = pd.read_csv(file, header=0)
    network_df.index = AA
    # Get chains from columns
    column_names = list(network_df.columns)
    chains = list(set([position[1] for position in column_names]))
    # Get positions without chain distinction from functional files
    positions = list(functional_data[protein].columns)
    average = pd.DataFrame(index=AA, columns=positions, dtype=np.float64)
    # Save data for position over chains in list, write average into df
    for position in positions:
        for aa in AA:
            values = []
            for chain in chains:
                check = position[0] + chain + position[1:]
                if check in network_df.columns:
                    values.append(network_df.at[aa, check])
            if values:
                average_value = sum(values)/len(values)
                average.at[aa, position] = average_value

    return average


def Standardize(protein, threshold, measure):
    """Return standardized values from network data. Make 0's into NaN.

    Given a protein, threshold and measure, data entry $x$  is standarized as:

    $$ x = \frac{x - \mu}{\sigma},$$

    where $\mu$ and $\sigma$ are the mean and standard deviation of the data.
    """

    network_df = ReadNetworkCSV(protein, threshold, measure)
    for position in network_df.columns:
        for aa in network_df.index:
            if position[0] == aa:
                network_df.at[aa, position] = np.nan
    data_array = network_df.to_numpy()
    data_mean = np.nanmean(network_df, dtype=np.float64)
    data_std = np.nanstd(network_df, dtype=np.float64)
    network_df = network_df.apply(lambda x:(x-data_mean)/data_std)

    return network_df


def GetPercentage(percentage, which, data, return_score=False):
    """Return set with top or bottom percentage of positions according to
    functional data.

    Parameters:
    ----------
    percentage (float): between 0 and 1, percentage of positions that we want.
    which (str): 'highest' (FRP), 'lowest' (FSP).
    data (dataframe): functional data to consider mean of
    return_score (bool): If True, return list of tuples with the mean
                         value per position.

    Returns:
    -------
    Set of positions FSP or FRP.
    """

    functional_mean = data.mean()
    positions = list(data.columns)
    pairs = [(functional_mean[pos], pos) for pos in positions]
    # Sort by increasing functional value
    pairs.sort(key = lambda x:x[0])
    if which == 'highest':
        pairs.reverse()
    n = int(len(positions) * percentage)
    if return_score:
        return [pair for pair in pairs[:n]]
    else:
        return set([pair[1] for pair in pairs[:n]])


def GetNetworkExtremes(protein, mincount, measure_cutoffs, thresh=9.0):
    """ Return the set of SSPs that pass measure sd cutoffs for at least
    mincount measures.

    Parameters
    ----------
    protein: str, PDB name of protein.
    mincount: int, minimal number of measures to consider.
    measure_cutoffs: list, cutoff vector.
    thresh: float, ditance threshold in Angstroms.
    """

    # In case thresh is given as an int
    thresh = float(thresh)
    network_extremes_list = []
    for i, measure in enumerate(measures[:3]):
        network_df = Standardize(protein, thresh, measure)
        if measure_cutoffs[i] > 0:
            extremes = network_df.columns[(network_df > measure_cutoffs[i]).any()].tolist()
        else:
            extremes = network_df.columns[(network_df < measure_cutoffs[i]).any()].tolist()
        network_extremes_list.extend(extremes)

    counter = Counter(network_extremes_list)
    positions = list(set(network_extremes_list))
    return set([pos for pos in positions if counter[pos] >= mincount])


def GetList(protein, mincount, measure_cutoffs, thresh=9.0, loss=True):
    """Get list with SSP positions, AA in three letter code.

    If loss==False, use complement for gain predictions. """

    pos = GetNetworkExtremes(protein, mincount, measure_cutoffs, thresh=thresh)
    if not loss:
        total_pos = functional_data[protein].columns
        complement = [i for i in total_pos if i not in pos]
        pos = complement
    positions = map(lambda x:pp.one_to_three(x[0])+x[1:], pos)

    return list(positions)


def ToPercentage(a,b):
    """Return percentage form of a/b, if b != 0. If given set or list, use len of.
    If string, return formatted percentage, else float."""
    x = a if type(a) == int or type(a) == float else len(a)
    y = b if type(b) == int or type(b) == float else len(b)

    if y == 0:
        return np.nan
    else:
        return round(100*x/y,1)
