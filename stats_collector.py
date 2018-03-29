"""collects statistics computed by feature_extractor.py and saves it to 1 .csv file - for future manipulation
"""
import os, math
from multiprocessing import Pool
from tqdm import trange
import numpy as np
import pandas as pd

def get_info(protein_path):
    info_file_path = os.path.join(protein_path, "info.txt")
    lines = []
    if os.path.exists(info_file_path):
        with open(info_file_path) as f:
            for line in f:
                line = line.strip()
                if len(line) > 0:
                    lines.append(line.replace("error", ";;") + "\n")
    else:
        lines = [";;\n"]
    return "".join([";".join([protein_path, line]) for line in lines])
    


def pool_processor(n):
    size = math.ceil(len(filenames)/total_processes)
    f = open("results_%s.csv"%n, 'w') 
    for index in trange(min(size, len(filenames)-(n-1)*size), position=n):
        data = get_info(filenames[(n-1)*size + index])
        f.write(data)
    f.close()


if __name__ == "__main__":
    total_processes = 10
    scPDBdir = "scPDB" # suppose that all data are located at this folder
    filenames = []
    with open("selected_subset.txt") as f:
        for line in f:
            filenames.append(os.path.join(scPDBdir, line.strip()))
    # next we open these folders and check if some data is present - in parallel way
    L = range(total_processes)
    with Pool(total_processes) as pool:
        pool.map(pool_processor, L)
    pool.join()
    print(1)
    df = []
    for i in L:
        print(i)
        filename = "results_%s.csv" % i
        df.append(pd.read_csv(filename, sep=";", header=None))
        os.remove(filename)
    df = pd.concat(df, ignore_index=True, axis=0)
    print(df.shape)
    df.to_csv("protein_subset_info.csv", sep=";", header=False, index=False)