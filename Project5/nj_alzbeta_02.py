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

#This is just Lena's code converted for Visual Studio
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
    minj = dist.columns[1]
    rimin = ris[0]
    rjmin = rjs[1]
    index_mini = 0
    index_minj = 1
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
                    rimin = ris[i]
                    rjmin = rjs[j]
                    index_mini = i
                    index_minj = j
    return(index_mini, index_minj, mini,minj,minN,rimin,rjmin)

def step2(mini,minj,dist):
    nodes = list(dist.index).copy()
    nodes.remove(mini)
    nodes.remove(minj)
    nodes.append('('+mini+','+minj+')')
    return nodes

def step3(index_mini, index_minj, mini, minj, dist, rimin, rjmin, node_length):
    nodes_length = list(dist.index).copy()
    nodes = node_length.copy()
    nodes.remove(node_length[index_mini])
    nodes.remove(node_length[index_minj])
    dik = round((dist.loc[mini,minj] + rimin - rjmin)/2,4)
    djk = round((dist.loc[mini,minj] + rjmin - rimin)/2,4)
    nodes.append('(' + node_length[index_mini] + ': '+ str(dik) + ','+ node_length[index_minj]+ ": "+ str(djk) + ')')
    return nodes

def step4(mini, minj, dist, nodes):
    k = nodes[-1]
    row_mini = dist.loc[mini,:]
    col_minj = dist[minj]
    for i in dist.index : 
        if i != mini and i!=minj : 
            dist.loc[i,k] = (dist.loc[mini,i]+dist.loc[minj,i]-dist.loc[mini,minj])/2
            dist.loc[k,i] = (dist.loc[mini,i]+dist.loc[minj,i]-dist.loc[mini,minj])/2
    dist.loc[k,k]=0
    dist = dist.drop(index = [mini,minj], columns = [mini,minj])
    return(dist)

def tree_construct_nj(dist):

    while(len(dist) > 3):

        index_mini, index_minj, mini,minj,minN,rimin,rjmin = step1(dist)
        nodes_1 = step2(mini, minj, dist)
        nodes_length = list(dist.index).copy()
        nodes_2 = step3(index_mini, index_minj, mini, minj, dist, rimin, rjmin, nodes_length)
        new_dist = step4(mini, minj, dist, nodes_1)
        print(new_dist)
        dist = new_dist


### main ###
dist = read_mtrx("C:\\Users\\alzbe\\Documents\\AU_Bioinfo_Masters\\Spring_2021\\AiB\\Projects\\Project_05\\example_slide4.phy")

tree_construct_nj(dist)