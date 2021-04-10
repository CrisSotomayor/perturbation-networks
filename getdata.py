import os
import Bio.PDB
import biographs as bg
import numpy as np
import networkx as nx
import csv


def GetData(path, prot, mutations, csv_path=None):
    """Get data from amino acid mutation perturbation networks as CSV files.

    Parameters:
        path (string): path where original and mutated pdb files are located
        prot (string): name of original pdb file
        mutations (dict): keys are strings representing positions to mutate
                        (amino acid, chain and index), each contains a list of
                        mutations performed (original aa, chain, index and
                        mutated aa). Mutated pdb files should be found in path,
                        as prot_mutation.pdb according to mutations in said list
        csv_path (string): default None, path where CSV files will be saved,
                        if None, a dir named "perturbation_network_data" will
                        be created in 'path'

    Returns:
        None
    """
    # Sorted list of one letter amino acids
    AA = list(Bio.PDB.Polypeptide.aa1)
    N = len(AA)  # Number of amino acids
    M = len(mutations.keys())  # Number of mutated positions
    cols = list(mutations.keys())  # List of mutated positions

    # Generate molecule of original pdb file
    original_prot = bg.Pmolecule(os.path.join(path, f"{prot}.pdb"))
    # The range of thresholds will define the networks
    thresholds = [round(i, 1) for i in np.linspace(3, 10, 71)]
    # Create dir to save resulting csv files if not specified
    if csv_path is None:
        csv_path = os.path.join(path, "perturbation_network_data")
        if not os.path.exists(csv_path):
            os.makedirs(csv_path)
    # Check if path and csv_path exist
    assert os.path.exists(path), f"Directory {path} doesn't exist."
    assert os.path.exists(csv_path), f"Directory {csv_path} doesn't exist."
    # Array to save data from original protein network
    original_data = np.zeros((4, len(thresholds)))

    # For each threshold we iterate over all mutations
    for i, threshold in enumerate(thresholds):
        nodes = np.zeros((N, M))
        edges = np.zeros((N, M))
        weight = np.zeros((N, M))
        distance = np.zeros((N, M))

        # Generate network for original graph with threshold
        original = original_prot.network(cutoff=threshold)
        original_matrix = nx.adjacency_matrix(original).toarray()
        # Saving data from original network
        original_data[0][i] = GetNodes(original)
        original_data[1][i] = GetEdges(original)
        original_data[2][i] = GetWeight(original)
        original_data[3][i] = GetDistance(original)

        for index, position in enumerate(mutations.keys()):
            for mutation in mutations[position]:
                # Generate network for current mutation
                current_path = os.path.join(path, f"{prot}_{mutation}.pdb")
                current_prot = bg.Pmolecule(current_path)
                current = current_prot.network(cutoff=threshold)

                # Obtain the absolute difference in terms of adjacency
                # matrices: the perturbation network.
                current_matrix = nx.adjacency_matrix(current).toarray()
                difference = np.abs(original_matrix - current_matrix)
                perturbation_network = nx.from_numpy_array(difference)

                # Remove isolates for accurate perturbation network node count
                perturbation_network.remove_nodes_from(
                    list(nx.isolates(perturbation_network)))

                # Corresponding row in array according to mutation
                assert mutation[-1] in AA, \
                    f"{mutation[-1]} not one of {Bio.PDB.Polypeptide.aa1}"
                aa_index = AA.index(mutation[-1])

                # Information obtained from perturbation network
                nodes[aa_index][index] = GetNodes(perturbation_network)
                edges[aa_index][index] = GetEdges(perturbation_network)
                weight[aa_index][index] = GetWeight(perturbation_network)
                distance[aa_index][index] = GetDistance(perturbation_network)

        # Save data arrays as csv files in csv_path
        WriteCSV(csv_path, nodes, cols, f"{prot}_{threshold}_nodes.csv")
        WriteCSV(csv_path, edges, cols, f"{prot}_{threshold}_edges.csv")
        WriteCSV(csv_path, weight, cols, f"{prot}_{threshold}_weight.csv")
        WriteCSV(csv_path, distance, cols, f"{prot}_{threshold}_distance.csv")

    # Save array from original data
    original_data = np.vstack(
        [thresholds, original_data])  # add thresholds
    original_data = np.transpose(original_data)  # to add names as header
    header = ['threshold', 'nodes', 'edges', 'weight', 'distance']
    # Write CSV of original data
    WriteCSV(csv_path, original_data, header, f"{prot}_original.csv")

    return None


def GetNodes(network):
    """ Return number of nodes in network. """
    return len(network.nodes())


def GetEdges(network):
    """ Return number of edges in network. """
    return len(network.edges())


def GetWeight(network):
    """ Return sum of weights for the edges in network. """
    return network.size(weight='weight')


def GetDistance(network):
    """ Return diameter of network, 0 if null graph, max of components if not connected. """
    if len(network.nodes()) == 0:  # If null graph, return 0
        return 0
    elif nx.is_connected(network) == False:  # If not connected, return maximun
        components = [network.subgraph(c).copy()
                      for c in nx.connected_components(network)]
        diameters = [nx.diameter(c) for c in components]
        return max(diameters)
    else:  # If connected, return diameter
        return nx.diameter(network)


def WriteCSV(csv_path, attribute, header, name):
    """ Write CSV file of attribute matrix with header in csv_path. """
    with open(os.path.join(csv_path, name), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(attribute)
    return None
