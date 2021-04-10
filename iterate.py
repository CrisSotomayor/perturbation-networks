""" Contains all code to run getmutations and getdata, iterating over individual
mutations instead of proteins """

import os
import sys
import csv
import Bio.PDB
import networkx as nx
import numpy as np
from collections import Counter
from Bio.PDB import MMCIFParser, PDBParser, Selection, NeighborSearch
import Bio.PDB.Polypeptide as pp

# Code from biographs bpdb.py
def pdb_model(structure_file, water=False, check_residues=True):
    """Return a biopython [1] model entity from a structure file.

    Parameters
    ----------
    structure_file: string
        Path to structure file
    water: boolean (default=False)
        True to take into account water molecules in the structure, False
        otherwise.
    check_residues: boolean (default=True)
        If True, removes all residues that are not standard amino acids.

    Notes
    -----
    1. http://biopython.org/wiki/Biopython

    """
    accepted_formats = ['cif', 'pdb', 'ent']
    parsers = [MMCIFParser, PDBParser, PDBParser]
    protein_name, file_format = structure_file.rsplit('.', 1)

    try:
        parser = parsers[accepted_formats.index(file_format)]
        parser = parser(QUIET=True)
    except ValueError:
        raise Exception("Accepted structure files are: {}".format(
            accepted_formats))

    structure = parser.get_structure(protein_name, structure_file)
    model = structure[0]

    if not water:
        for chain in model.get_chains():
            for residue in list(chain):
                hetero_flag = residue.id[0].strip()
                # Empty strings evaluate to False.  Therefore hetero_flag
                # returns False if the residue is not a water molecule.
                if hetero_flag:
                    chain.detach_child(residue.id)
            if not list(chain):
                model.detach_child(chain.id)

    if check_residues:
        for chain in model.get_chains():
            for residue in list(chain):
                # FoldX removes hetatom flags, remove leftover ions and non standard
                # amino acids
                if not pp.is_aa(residue, standard=True):
                    chain.detach_child(residue.id)
                # Fix for issue where single N atom is read as aa
                elif len(residue.get_unpacked_list()) == 1:
                    chain.detach_child(residue.id)
            if not list(chain):
                model.detach_child(chain.id)

    return model

def label_residue(residue):
    """ Return a string of the label of the biopython [1] residue object.

    The label of the residue is the following:
        Chain + Position

    Parameters
    ----------
    residue: Bio.PDB.Residue.Residue
        The residue to be labeled.

    Notes
    -----
    1. http://biopython.org/wiki/Biopython

    """
    position = str(residue.id[1])
    chain = residue.parent.id

    return chain + position

def residue_adjacency(model, cutoff=5, weight=True):
    """Return residue adjacency dictionary defined by cutoff distance.

    Parameters
    ----------
    model: Bio.PDB.Model
        Model created with the atomic coordinates of the protein file.

    cutoff: int or float
        Distance cutoff defining links between atoms.  Two atoms are adjacent
        if their distance is less than the given cutoff.

    See Also
    --------
    pdb_model

    """
    atoms = Selection.unfold_entities(model, 'A')

    neighbor_search = NeighborSearch(atoms)
    atomic_adjacency = {}

    for atom in atoms:
        _res = label_residue(atom.get_parent())
        adjacent_atoms = []
        for adj_atom in neighbor_search.search(atom.coord, cutoff):
            _adj_res = label_residue(adj_atom.parent)
            # Adjacent atoms must be in different residues
            if _adj_res != _res:
                adjacent_atoms.append(adj_atom)
        atomic_adjacency[atom] = adjacent_atoms

    adjacency = {}

    # Create residue adjacency dictionary with string format, see
    # label_residue.
    for atom, neighbors in atomic_adjacency.items():
        residue = label_residue(atom.get_parent())
        adjacency.setdefault(residue, [])

        # Only different residues are connected by an edge (No loops).
        not_in_residue = []
        for neighbor in neighbors:
            neighbor_parent = label_residue(neighbor.get_parent())
            if neighbor_parent is not residue:
                not_in_residue.append(neighbor_parent)

        adjacency[residue].extend(not_in_residue)

    if not weight:

        return adjacency

    # Make new dictionary mapping each residue to its neighbors taking
    # into account the weight.
    weighted_adjacency = {}
    for residue in adjacency:
        counter = Counter(adjacency[residue])
        weighted_adjacency[residue] = {
            neighbor: {'weight': counter[neighbor]}
            for neighbor in counter}

    return weighted_adjacency

# Code from biographs bgraph.py
def network(model, cutoff=5, weight=True):
    """Return the interaction network of a protein structure.

    The interaction network is defined by a distance cutoff.

    Parameters
    ----------
    model: Bio.PDB.model
        The protein structure.
    cutoff: float
        The distance cutoff defining an interaction between two nodes.
    weight: boolean
        True if atomic interactions are to be considered.
    """

    adjacency_dictionary = residue_adjacency(model, cutoff=cutoff,
                                            weight=weight)

    return nx.Graph(adjacency_dictionary)

# Code from biographs Pmolecule.py
class Pmolecule(object):
    """Create a Pmolecule object.

    The Pmolecule calls a number of methods for the analysis of protein
    structure. This includes the contruction of the interaction network of the
    protein.

    Parameters
    ----------
    structure_file : str
        The path to the structure file of the targeted protein. Three
        structure-file formats are supported: `pdb', `cif', and `ent'.
    water : boolean, (default: False)
        If False, water molecules are removed.
    check_residues : boolean, (default: True)
        If True, removes all residues that are not standard amino acids.

    Attributes
    ----------
    model: Bio.PDB.model
        The structural model of the structure. See www.biopython.org.
    network: networkx:Graph
        The amino acid network of the protein based on a distance
        cutoff (default=5 angs). See :Pmolecule:network: for more info.
    path_to_file: str
        The path to the structural file used to instantiate the class.
    """

    def __init__(self, structure_file, water=False, check_residues=True):

        self.model = pdb_model(structure_file, water=water, check_residues=check_residues)
        self.path_to_file = structure_file

    def __len__(self):
        """Returns number of residues in the molecule"""

        return len(list(self.model.get_residues()))

    def residue(self, node):
        """Return :Bio:PDB:Residue: from node of network

        node: str, node of the form :chain:position:
        """

        if type(node) is not str:
            raise Exception(
                "{} should be a string of the form "
                ":chain:position:".format(node))
        try:
            res = self.model[node[0]][int(node[1:])]
        except KeyError:
            raise Exception("node {} is not a residue in molecule".format(
                node))

        chain = res.parent
        pos = res.id[1]

        return self.model[chain.id][pos]

    def network(self, cutoff=5, weight=True):
        """Return the interaction network of a protein structure.

        The interaction network is defined by a distance cutoff.

        Parameters
        ----------
        model: Bio.PDB.model
            The protein structure.
        cutoff: float
            The distance cutoff defining an interaction between two nodes.
        weight: boolean
            True if atomic interactions are to be considered.
        """

        return network(self.model, cutoff=cutoff, weight=weight)

# Get list of mutations to iterate over
def MutationsList(path, protein):
    """Get list with all possible mutations of protein, keeps mutations that
    exist as pdbs in path.

    Parameters:
        path (string): path where protein pdb file is located, and mutated pdbs
                       if pdbs is True
        protein (string): name of original pdb file

    Returns:
        mutations: list with :aa:chain:position:mutated_aa for all mutations
        positions = list with :aa:chain:position for mutated positions
    """
    # Sorted list of one letter amino acids
    AA = list(pp.aa1)
    # Generate model of original pdb file
    model = Pmolecule(os.path.join(path, f"{protein}.pdb")).model

    # List to store mutations
    mutations = []
    # List to store mutated positions
    positions = []

    for chain in model.get_chains():
        for residue in chain:
            if pp.is_aa(residue):
                code = pp.three_to_one(residue.get_resname())
                chain = residue.parent.id
                position = str(residue.id[1])
                prefix = code+chain+position
                if prefix[0] != AA[0]:
                    first = AA[0]
                else:
                    first = AA[1]
                first_mutation = os.path.join(path, f"{protein}_{prefix+first}.pdb")
                # Assume if first mutation exists, all mutations do
                # Excludes mutating residue for itself
                if os.path.exists(first_mutation):
                    mutations.extend([prefix+aa for aa in AA if aa!=prefix[0]])
                    positions.append(prefix)

    return mutations, positions

# Get adjacency matrix from residue_adjacency
def AdjacencyMatrix(model, cutoff=5):
    """Get adjacency matrix from residue_adjacency dictionary.

    Parameters:
        model: model returned from pdb_model

    Returns:
        adjacency_matrix: ndarray, adjacency matrix for model with cutoff
    """
    # List with all residues in model, chain:index
    positions = [label_residue(residue) for residue in model.get_residues()]

    N = len(positions)
    adjacency_matrix = np.zeros((N,N))

    weighted_adjacency = residue_adjacency(model, cutoff=cutoff)

    for residue in weighted_adjacency:
        for neighbor in weighted_adjacency[residue]:
            i = positions.index(residue)
            j = positions.index(neighbor)
            adjacency_matrix[i][j] = weighted_adjacency[residue][neighbor]['weight']

    return adjacency_matrix

# Functions to get data, from getdata.py
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
    if not os.path.exists(csv_path):
        os.makedirs(csv_path)
    with open(os.path.join(csv_path, name), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        writer.writerows(attribute)
    return None

# New get data function, modified from GetData in getdata.py
def GetData(protein, mutation, path, threshold):
    """Get nodes, edges, weight and distance from mutation.

    Parameters:
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

    # Generate network for current mutation
    current_path = os.path.join(path, f"{protein}_{mutation}.pdb")
    current_prot = Pmolecule(current_path)
    current_model = current_prot.model

    # Obtain the absolute difference in terms of adjacency
    # matrices: the perturbation network.
    current_matrix = AdjacencyMatrix(current_model, cuttoff=threshold)
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
    nodes[aa_index][index] = GetNodes(perturbation_network)
    edges[aa_index][index] = GetEdges(perturbation_network)
    weight[aa_index][index] = GetWeight(perturbation_network)
    distance[aa_index][index] = GetDistance(perturbation_network)

    return
