from time import sleep
from scipy.spatial.distance import euclidean
import pandas as pd
from multiprocessing import Pool
import numpy as np
import os, math
from tqdm import trange
import h5py

from pyvdwsurface import hmsurface
from biopandas.mol2 import PandasMol2


def rotate_relative_to_ligand(coords, ligand_coords):
    mean_ligand_coords = np.mean(ligand_coords, axis=0)
    mean_coords = np.mean(coords, axis=0)
    ligand_centered = ligand_coords - mean_ligand_coords
    u, s, vh = np.linalg.svd(ligand_centered, full_matrices=False)
    ligand_rotated = np.dot(ligand_centered, vh.T)
    coords_rotated = np.dot(coords - mean_ligand_coords, vh.T)
    return coords_rotated, ligand_rotated


def get_subst_description(mol2_lines):
    started = False
    for i, x in enumerate(mol2_lines):

        if x.startswith("@<TRIPOS>SUBSTRUCTURE"):
            started = True
            start_index = i+1
        elif started:
            if x.startswith("@<TRIPOS>"):
                break
            else:
                end_index = i
    return mol2_lines[start_index: end_index+1]


def get_substructure_df(mol2_lines):
    lst1 = get_subst_description(mol2_lines)
    column_names = 'subst_id subst_name number type chain_id chain_name aa_name'.split()
    column_types = (int, str, int, str, int, str, str)

    import pandas as pd
    df = pd.DataFrame([lst.split()[:len(column_names)] for lst in lst1],
                              columns=column_names)
    for i in range(df.shape[1]):
        df[column_names[i]] = df[column_names[i]].astype(column_types[i])
    return df

def process_protein(path, density=0.5):
    site_file = os.path.join(path, "site.mol2")
    ligand_file = os.path.join(path, "ligand.mol2")
    site = PandasMol2().read_mol2(site_file)
    
    ligand = PandasMol2().read_mol2(ligand_file)
    df = get_substructure_df(site.mol2_text.split('\n'))
    subst_names = df[df["type"]=="RESIDUE"]["subst_name"].values
    filtered_site = site.df[site.df["subst_name"].isin(subst_names)]
    atom_names = filtered_site["atom_type"].apply(lambda x: x.split(".")[0].encode())
    coords = np.ascontiguousarray(filtered_site["x y z".split()].values, dtype=np.float64)
    ligand_coords = np.ascontiguousarray(ligand.df["x y z".split()].values, dtype=np.float64)
    c, l = rotate_relative_to_ligand(coords, ligand_coords)
    surface_points = hmsurface(c, atom_names, l, density=density)
    return surface_points, l


def save_protein_info(basedir, protein, savedir, density=0.5):
    datapath = os.path.join(basedir, protein)
    savepath = os.path.join(savedir, protein+".hdf5")
    with h5py.File(savepath, 'w') as f:
        try:
            points, ligand = process_protein(datapath, density=density)
        except:
            with open("problems.txt", 'a') as f:
                f.write(protein+" - exception found during processing\n")
            return
        f.create_dataset("points", data=points, compression="gzip")
        f.create_dataset("ligand", data=points, compression="gzip")
        f.creade_dataset("density", data=density)


def get_proteins_list(annotations_file="scPDB_Results.tsv"):
    import pandas as pd
    df = pd.read_table(annotations_file)
    simple_df = df[df["ChainPercentageInSite"].apply(lambda x: x.find("//")<0)]
    return simple_df.apply(lambda row:row["PDB_ID"] + "_"+str(row["Site_Number"]), axis=1).values


if __name__ == "__main__":
    scPDBdir = "/home/malygina/data/scPDB" # suppose that all data are located at this folder
    filenames = get_proteins_list()
    path2save = "/home/malygina/data/generated/0_5"
    density = 0.5
    if not os.path.exists(path2save):
        os.mkdir(path2save)
        
    total_processes = 10
    L = range(total_processes)

    def pool_processor(n):
        size = math.ceil(len(filenames)/total_processes)
        for index in trange(min(size, len(filenames)-(n-1)*size), position=n):
            save_protein_info(scPDBdir, filenames[(n-1)*size + index], path2save)
 
    with Pool(len(L)) as pool:
        pool.map(pool_processor, L)