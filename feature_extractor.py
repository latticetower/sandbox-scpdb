from htmd.ui import *
from htmd.molecule.voxeldescriptors import getVoxelDescriptors
from scipy.spatial.distance import euclidean


def get_windows_data(data, ndims=(16, 16, 16)):
    """helper iterator method: gets data, slices it to windows defined by ndims, returns its coordinates 
    """
    if isinstance(ndims, tuple):
        (nx, ny, nz) = ndims
    else:
        nx = ny = nz = ndims
    (maxx, maxy, maxz) = data.shape[:3]
    # the following might be ineffective, but simple
    for i in range(maxx - nx):
        for j in range(maxy - ny):
            for k in range(maxz - nz):
                yield i, j, k, data[i:i+nx, j:j+ny, k:k+nz, :]

        
def process_folder_balanced(protein_path, ndims=16, cutoff=4):
    """we return subset with positive and neg values, balanced by downsampling
    """
    #total_features = []
    binding_site = Molecule(os.path.join(protein_path, "site.mol2")) 
    binding_site_center = np.mean(binding_site.get('coords'), axis=0)
    #print(binding_site_center)
    protein_molecule = Molecule(os.path.join(protein_path, "protein.mol2"))
    features, centers, N = getVoxelDescriptors(protein_molecule, buffer=8)
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
    neg_downsampled = np.random.choice(neg_values, pos_values.shape[0])
    coords = np.concatenate([coords[pos_values], coords[neg_downsampled]])
    y = np.concatenate([y[pos_values], y[neg_downsampled]])
    # next we permute subset
    ids = np.random.permutation(y.shape[0])
    for index in ids:
        (i, j, k) = coords[index]
        yield features[i: i+ndims, j: j+ndims, k: k+ndims, :], y[index]

    
if __name__ == "__main__":
    scPDBdir = "scPDB" # suppose that all data are located at this folder
    ## full usage example:
    # for protein_folder_name in os.listdir(scPDBdir):
    #    protein_folder_path = os.path.join(scPDBdir, protein_folder_name)
    #    for features, y in process_folder(protein_folder_path, cutoff=4):
    #        print("do something") # here we have a list of preprocessed features of size (16, 16, 16, 8) and y labels
    #
    # for features, y in process_folder("scPDB/4f9g_1", cutoff=4):
    #    print(y) # or iterate over all
    for features, y in process_folder_balanced("scPDB/4f9g_1", cutoff=4):
        print(features.shape, y) # or do something
    pass