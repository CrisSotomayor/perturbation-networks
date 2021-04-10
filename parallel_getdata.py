import os
import multiprocessing as mp
import iterate as it
import numpy as np
import networkx as nx
import Bio.PDB
from ctypes import c_char_p
from multiprocessing import sharedctypes



def GetData(tuple):
    """Get nodes, edges, weight and distance from mutation.

    Parameters:
        tuple containing:
            protein (str): name of original protein pdb file
            mutation (str): :aa:chain:position:mutated_aa
            path (str): path were pdbs are found, assumes mutated file is named
                        'protein_mutation.pdb'
            threshold (float): threshold to generate network

    Uses:
        original_matrix = matrix from original protein
        nodes, edges, weighst, distance = np arrays to write data in
        positions = list of positions, to get corresponding index to write in


    Returns: None
    """
    protein = tuple[0]
    mutation = tuple[1]
    path = tuple[2]
    threshold = tuple[3]

    # Generate network for current mutation
    current_path = os.path.join(path, f"{protein}_{mutation}.pdb")
    current_prot = it.Pmolecule(current_path)
    current_model = current_prot.model

    # Obtain the absolute difference in terms of adjacency
    # matrices: the perturbation network.
    current_matrix = it.AdjacencyMatrix(current_model, cutoff=threshold)
    original_matrix_np = np.ctypeslib.as_array(original_matrix)
    difference = np.abs(original_matrix_np - current_matrix)
    perturbation_network = nx.from_numpy_array(difference)

    # Remove isolates for accurate perturbation network node count
    perturbation_network.remove_nodes_from(
        list(nx.isolates(perturbation_network)))

    # Corresponding row in array according to mutation
    assert mutation[-1] in AA, \
        f"{mutation[-1]} not one of {Bio.PDB.Polypeptide.aa1}"
    aa_index = AA.index(mutation[-1])

    # Corresponding column in array according to position
    assert mutation[:-1] in positions, \
        f"{mutation[:-1]} not one of amino acids in the protein"
    index = positions.index(mutation[:-1])

    # Information obtained from perturbation network
    nodes[aa_index][index] = it.GetNodes(perturbation_network)
    edges[aa_index][index] = it.GetEdges(perturbation_network)
    weight[aa_index][index] = it.GetWeight(perturbation_network)
    distance[aa_index][index] = it.GetDistance(perturbation_network)

    return

if __name__ == '__main__':

    thresholds = [round(i, 1) for i in np.linspace(3, 10, 71)]
    proteins = ["4bz3", "1d5r"]

    home_path = os.path.join(os.getenv("HOME"), "perturbation_networks")

    with mp.Manager() as manager:
        for protein in proteins:
            # Path where all necessary pdbs are stored
            path = os.path.join(home_path, f"{protein}")
            # Path were resultings data csvs will be stored
            csv_path = os.path.join(home_path, "data")
            # Dir in path for protein
            prot_csv_path = os.path.join(csv_path, f"{protein}")
            # Create paths if they don't exist
            if not os.path.exists(csv_path):
                os.makedirs(csv_path)
            if not os.path.exists(prot_csv_path):
                os.makedirs(prot_csv_path)
            # Original data matrix
            original_data = np.zeros((4, len(thresholds)))

            AA = manager.list(list(Bio.PDB.Polypeptide.aa1)) # List of amino acids
            mutations_positions = it.MutationsList(path, protein)
            mutations = manager.list(mutations_positions[0])
            positions = manager.list(mutations_positions[1])
            N = len(AA) # Number of amino acids
            M = len(positions) # Number of mutated positions

            for i, threshold in enumerate(thresholds):
                # Create arrays to store data
                # https://jonasteuwen.github.io/numpy/python/multiprocessing/2017/01/07/multiprocessing-numpy-array.html
                nodes1 = np.ctypeslib.as_ctypes(np.zeros((N, M)))
                edges1 = np.ctypeslib.as_ctypes(np.zeros((N, M)))
                weight1 = np.ctypeslib.as_ctypes(np.zeros((N, M)))
                distance1 = np.ctypeslib.as_ctypes(np.zeros((N, M)))
                # Shared arrays
                nodes =  sharedctypes.RawArray(nodes1._type_, nodes1)
                edges =  sharedctypes.RawArray(edges1._type_, edges1)
                weight =  sharedctypes.RawArray(weight1._type_, weight1)
                distance =  sharedctypes.RawArray(distance1._type_, distance1)

                # Generate molecule of original pdb file
                original_prot = it.Pmolecule(os.path.join(path, f"{protein}.pdb"))
                original_model = original_prot.model
                original = original_prot.network(cutoff=threshold)
                matrix = np.ctypeslib.as_ctypes(it.AdjacencyMatrix(original_model, cutoff=threshold))
                # Shared data for original matrix with threshold
                original_matrix = sharedctypes.RawArray(matrix._type_, matrix)

                # Save original data
                original_data[0][i] = it.GetNodes(original)
                original_data[1][i] = it.GetEdges(original)
                original_data[2][i] = it.GetWeight(original)
                original_data[3][i] = it.GetDistance(original)

                # Get list of tuples to map in Pool
                tuples = [(protein, mutation, path, threshold) for mutation in mutations]

                with mp.Pool(30) as pool:
                    pool.map(GetData, tuples)

                # Arrays should be full by now, we save as csvs in folder in csv_path
                # positions will be used as header
                it.WriteCSV(prot_csv_path, nodes, positions, f"{protein}_{threshold}_nodes.csv")
                it.WriteCSV(prot_csv_path, edges, positions, f"{protein}_{threshold}_edges.csv")
                it.WriteCSV(prot_csv_path, weight, positions, f"{protein}_{threshold}_weight.csv")
                it.WriteCSV(prot_csv_path, distance, positions, f"{protein}_{threshold}_distance.csv")

            # Save array from original data
            original_data = np.vstack(
                [thresholds, original_data])  # add thresholds
            original_data = np.transpose(original_data)  # to add names as header
            header = ['threshold', 'nodes', 'edges', 'weight', 'distance']
            # Write CSV of original data
            it.WriteCSV(prot_csv_path, original_data, header, f"{protein}_original.csv")
