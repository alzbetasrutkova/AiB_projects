import sys
import numpy as np
import pandas as pd

def read_mtrx(filename):
    f = open(filename,'r')
    num_species = int(f.readline())

    dict_mat = {}
    mat = np.zeros((num_species,num_species))
    
    for i in range(0,num_species):
        line = f.readline()
        nums_in_line = line.split()
        dict_mat[nums_in_line[0]] = i
        for j in range(1,num_species +1):
            mat[i,j-1] = nums_in_line[j]
    f.close()

    species = list(dict_mat.keys())
    dict_dist = {}
    for i in range(len(species)) :
        dict_dist[species[i]] = mat[i]
    dict_dist

    dist=pd.DataFrame.from_dict(dict_dist, orient='index', columns=species)
    return dist

def step1(dist):
    S = len(dist.columns)
    N = pd.DataFrame(columns=dist.columns, index = dist.index)
    ris = []
    for i in dist.index:
        ri = 0
        for j in dist.columns:
            ri +=dist.loc[i,j]
        ris.append(ri/(S-2))
    
    rjs = []
    for j in dist.columns:
        rj = 0
        for i in dist.index:
            rj +=dist.loc[i,j]
        rjs.append(rj/(S-2))
        
    minN = sys.maxsize
    mini = dist.index[0]
    minj = dist.columns[0]
    for i in range(len(dist.index)):
        for j in range(len(dist.columns)):
            ci = dist.index[i]
            cj = dist.columns[j]
            nij = dist.loc[ci,cj]-ris[i]-rjs[j]
            N.loc[ci,cj]=nij
            if ci != cj:
                if nij<minN : 
                    minN = nij
                    mini = ci
                    minj = cj
    return(mini,minj,minN)

### main ###
dist = read_mtrx("C:\\Users\\alzbe\\Documents\\AU_Bioinfo_Masters\\Spring_2021\\AiB\\Projects\\Project_05\\example_slide4.phy")
print(dist)
print(step1(dist))