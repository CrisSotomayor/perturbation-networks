import os
import subprocess
import shutil
import Bio.PDB
import Bio.PDB.Polypeptide as pp
#import biographs as bg
import iterate as it
import multiprocessing as mp

def MutationsDict(file, positions=None):
    """Get dictionary with lists of mutations per position in protein, ignore
    positions without residue in pdb file.

    Parameters:
        file (string): pdb file to get mutations from
        positions: list of tuples of the form (chain, first, last) for positions
                   to mutate for all other aminoacids. If None, mutates all
                   positions in all chains

    Returns:
        dict with keys :aa:chain:position, each containing lists with
        :aa:chain:position:mutated_aa for all mutations

    """

    # Sorted list of one letter amino acids
    AA = list(Bio.PDB.Polypeptide.aa1)
    # Generate model of original pdb file
    model = it.Pmolecule(file).model
    # Dict to store mutations
    mutations = dict()
    if positions:
        for chain_id, first, last in positions:
            # Get chain corresponding to chain_id given
            chain = next(chain for chain in model.get_chains() if chain.id == chain_id)
            for residue in chain:
                if pp.is_aa(residue):
                    code = pp.three_to_one(residue.get_resname())
                    position = residue.id[1]
                    prefix = code+chain_id+str(position)
                    # Only save positions between first and last
                    if position in range(first, last +1):
                            mutations[prefix] = [prefix+aa for aa in AA if aa!=code]
    else:
        for chain in model.get_chains():
            for residue in chain:
                if pp.is_aa(residue):
                    code = pp.three_to_one(residue.get_resname())
                    position = residue.id[1]
                    chain_id = chain.id
                    prefix = code+chain_id+str(position)
                    mutations[prefix] = [prefix+aa for aa in AA if aa!=code]
    return mutations

def GetMutations(path, protein, mutations, foldx_path, out_path=None, use_mp=None):
    """Get mutated pdb files using FoldX.

    Parameters:
        path (string): path where original pdb file is stored, if out_path is None,
                       all resulting files will be stored here
        protein (string): name of original pdb file, without .pdb suffix
        mutations (dict): contains mutations to be made, as returned by MutationsDict
        foldx_path (string): path where FoldX software is stored
        out_path (string): path where resulting files will be stored. Default None
        use_mp (int): number of pools to run if multiprocessing, if None, multiprocessing
                  is not used

    Returns:
        None, only output files

    """
    positions = mutations.keys()
    assert os.path.exists(path), f"{path} does not exist."
    if not out_path:
        out_path = path
    else:
        if not os.path.exists(out_path):
            os.makedirs(out_path)

    # Mutations will be made according to :aa:chain:position items in positions,
    # generating an individual list and config file for each, which are necessary
    # to run FoldX BuildModel
    if use_mp:
        for position in positions:
            # Since output FoldX files are named protein_i, we need to store
            # them separately so they don't overwrite each other
            current_path = os.path.join(out_path, position)
            os.makedirs(current_path, exist_ok=True)
            IndList(current_path, position, mutations[position])
            ConfigFile(path, current_path, protein, position)
        with mp.Manager() as manager:
            mutations = manager.dict(mutations)
            inputs = [(os.path.join(out_path, position), protein, position,
                        mutations[position], foldx_path) for position in positions]
            with mp.Pool(use_mp) as pool:
                pool.starmap(Mutate, inputs)
    else:
        for position in positions:
            IndList(out_path, position, mutations[position])
            ConfigFile(path, out_path, protein, position)
            Mutate(out_path, protein, position, mutations[position], foldx_path, use_mp=False)
    return

def Mutate(path, protein, prefix, mutations_list, foldx_path, use_mp=True):
    """Call FoldX BuildModel through config_file, rename resulting .pdb files,
    delete files not needed."""
    config_file = os.path.join(path, f"config_{protein}_{prefix}.cfg")
    mutate = f"{foldx_path} -f {config_file}"
    subprocess.check_call(mutate, shell = True)
    for i, mutation in enumerate(mutations_list):
        if use_mp:
            # Rename and move, assume we're in path/position, move pdb files to path
            parent_path = parent_path = os.path.split(path)[0]
            source = os.path.join(path, f"{protein}_{i+1}.pdb")
            dest = os.path.join(parent_path, f"{protein}_{mutation}.pdb")
            os.rename(source, dest)
        else:
            # Rename file according to mutation
            source = os.path.join(path, f"{protein}_{i+1}.pdb")
            dest = os.path.join(path, f"{protein}_{mutation}.pdb")
            os.rename(source, dest)
        # Delete corresponding WT file
        WT_file = os.path.join(path, f"WT_{protein}_{i+1}.pdb")
        os.remove(WT_file)
    # Delete all other resulting FoldX files, config and ind_list
    file_names = [f"Average_{protein}_{prefix}_{protein}.fxout",
                  f"config_{protein}_{prefix}.cfg",
                  f"Dif_{protein}_{prefix}_{protein}.fxout",
                  f"individual_list_{prefix}.txt",
                  f"PdbList_{protein}_{prefix}_{protein}.fxout",
                  f"Raw_{protein}_{prefix}_{protein}.fxout"]
    to_delete = [os.path.join(path, file) for file in file_names]
    for file in to_delete:
        os.remove(file)
    if use_mp:
        # Delete leftover empty dir
        os.rmdir(path)
    return

def IndList(path, prefix, mutation_list):
    """Write individual_list txt file."""
    file = os.path.join(path, f"individual_list_{prefix}.txt")
    indlist = open(file, "w+")
    for mutation in mutation_list:
        indlist.write(f"{mutation};\n")
    indlist.close()
    return

def ConfigFile(protein_path, out_path, protein, prefix):
    """Write config file."""
    file = os.path.join(out_path, f"config_{protein}_{prefix}.cfg")
    config = open(file, "w+")
    config.write("command=BuildModel\n")
    config.write(f"pdb={protein}.pdb\n")
    config.write(f"pdb-dir = {protein_path}\n")
    ind_list = os.path.join(out_path, f"individual_list_{prefix}.txt")
    config.write(f"mutant-file={ind_list}\n")
    config.write(f"output-dir={out_path}\n")
    config.write(f"output-file={protein}_{prefix}\n")
    config.close()
    return

if __name__ == '__main__':
    protein = '1d5r'
    home_path = os.path.join(os.getenv("HOME"), "perturbation_networks")
    # Path where original pdb and foldx software are stored
    path = os.path.join(home_path, protein)
    foldx_path = os.path.join(home_path, "foldx/foldx")
    pdb_file = os.path.join(path, f"{protein}.pdb")
    mutations = MutationsDict(pdb_file)
    GetMutations(path, protein, mutations, foldx_path, out_path=path, use_mp=30)
