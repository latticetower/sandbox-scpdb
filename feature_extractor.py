from time import sleep
from htmd.ui import *
from htmd.molecule.voxeldescriptors import getVoxelDescriptors
from scipy.spatial.distance import euclidean
import pandas as pd
from multiprocessing import Pool
import numpy as np
import os, math
from tqdm import trange


def get_windows_data(data, ndims=(16, 16, 16)):
    """helper iterator method: gets data, slices it to windows defined by ndims, returns its coordinates 
    """
    if isinstance(ndims, tuple):
        (nx, ny, nz) = ndims
    else:
        nx = ny = nz = ndims
    (maxx, maxy, maxz) = data.shape[:3]
    # the following might be ineffective, but simple
    for i in range(0, maxx - nx, 4):
        for j in range(0, maxy - ny, 4):
            for k in range(0, maxz - nz, 4):
                yield i, j, k, data[i:i+nx, j:j+ny, k:k+nz, :]

        
def process_folder_balanced(protein_path, ndims=16, cutoff=4):
    """we return subset with positive and neg values, balanced by downsampling
    """
    info_file_path = os.path.join(protein_path, "info.txt")
    for x in ["site.mol2", "protein.mol2"]:
        if not os.path.exists(os.path.join(protein_path, x)):
            return
    #total_features = []
    binding_site = Molecule(os.path.join(protein_path, "site.mol2")) 
    binding_site_center = np.mean(binding_site.get('coords'), axis=0)
    #print(binding_site_center)
    protein_molecule = Molecule(os.path.join(protein_path, "protein.mol2"))
    try:
        features, centers, N = getVoxelDescriptors(protein_molecule, buffer=8)
    except:
        with open(info_file_path, 'w')as f:
            f.write("error")
        return
    if np.any(N < ndims):
        print("protein is smaller, skip")
        return
    N = tuple(N) + (-1,)
    features = features.reshape(*N)
    centers = centers.reshape(*N)
    coords = [] # coords of the starting point for (16, 16, 16) block
    y = []
    for i, j, k, window_data in get_windows_data(centers, ndims=ndims):
        voxel_center = np.mean(window_data.reshape((-1, window_data.shape[-1])), axis=0)
        coords.append((i, j, k))        
        y.append(euclidean(voxel_center, binding_site_center) < cutoff)
    coords = np.asarray(coords)
    y = np.asarray(y)
    pos_values = np.where(y)[0]
    neg_values = np.where(y==0)[0]
    with open(info_file_path, 'w') as f:
        f.write("%s;%s;%s"% (len(y), len(pos_values), len(neg_values)))
        
    neg_downsampled = np.random.choice(neg_values, pos_values.shape[0])
    coords = np.concatenate([coords[pos_values], coords[neg_downsampled]])
    y = np.concatenate([y[pos_values], y[neg_downsampled]])
    # next we permute subset
    ids = np.random.permutation(y.shape[0])
    for index in ids:
        (i, j, k) = coords[index]
        yield features[i: i+ndims, j: j+ndims, k: k+ndims, :], y[index]


def process_protein_data(protein_path, recompute=False):
    ndims=16
    csv_filename = os.path.join(protein_path, "features.csv")
    if not recompute and os.path.exists(csv_filename):
        return
    data = pd.DataFrame(columns=["y"] + list(range(0, ndims*ndims*ndims*8)))
    for i, (features, y) in enumerate(process_folder_balanced(protein_path, ndims=ndims, cutoff=4)):
        data.loc[i, "y"] = y
        data.loc[i, 1:] = features.reshape(-1)
    data.to_csv(csv_filename, sep=";", index=False, header=False)
   
    
if __name__ == "__main__":
    scPDBdir = "scPDB" # suppose that all data are located at this folder
    filenames = []
    with open("selected_subset.txt") as f:
        for line in f:
            filenames.append(os.path.join(scPDBdir, line.strip()))
    #filenames = list(map(lambda x: os.path.join(scPDBdir, x), os.listdir(scPDBdir)))
    #filenames = filenames[:]
    total_processes = 10
    L = range(total_processes)

    def pool_processor(n):
        size = math.ceil(len(filenames)/total_processes)
        for index in trange(min(size, len(filenames)-(n-1)*size), position=n):
            process_protein_data(filenames[(n-1)*size + index])
 
    with Pool(len(L)) as pool:
        pool.map(pool_processor, L)
