## NJ
import numpy as np
import pandas as pd
import sys

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
    return dict_mat, mat 

def step1(dist):
    S = len(dist.columns)
    N = pd.DataFrame(columns=dist.columns, index = dist.index)
    rs = []
    for i in dist.index:
        ri = 0
        for j in dist.columns:
            ri +=dist.loc[i,j]
        rs.append(ri/(S-2))
            
    minN = sys.maxsize
    mini = dist.index[0]
    minj = dist.columns[1]
    rimin = rs[0]
    rjmin = rs[1]
    index_mini = 0
    index_minj = 1
    for i in range(len(dist.index)):
        for j in range(len(dist.columns)):
            ci = dist.index[i]
            cj = dist.columns[j]
            nij = dist.loc[ci,cj]-rs[i]-rs[j]
            N.loc[ci,cj]=nij
            if ci != cj:
                if nij<minN : 
                    minN = nij
                    mini = ci
                    minj = cj
                    rimin = rs[i]
                    rjmin = rs[j]
                    index_mini = i
                    index_minj = j
    return(index_mini, index_minj, mini,minj,minN,rimin,rjmin)

def NJ(dist_mat):
    # create the distance dataframe
    species_dict, dist = read_mtrx(dist_mat)
    species = list(species_dict.keys())
    dict_dist = {}
    for i in range(len(species)) :
        dict_dist[species[i]] = dist[i]
    dist=pd.DataFrame.from_dict(dict_dist, orient='index', columns=species) 
    
    # initializing variables
    nodes_length = list(dist.index).copy()
    nodes = list(dist.index).copy()
    
    while len(dist.index)>3 : 
        
        #step 1 
        index_mini, index_minj,mini,minj, minN,rimin, rjmin = step1(dist)
        
        # step 2 
        nodes = list(dist.index).copy()
        nodes.remove(mini)
        nodes.remove(minj)
        nodes.append('('+mini+','+minj+')')
        
        #step 3 
        nodes2 = nodes_length.copy()
        nodes_length.remove(nodes2[index_mini])
        nodes_length.remove(nodes2[index_minj])
        dik = round((dist.loc[mini,minj] + rimin - rjmin)/2,4)
        djk = round((dist.loc[mini,minj] + rjmin - rimin)/2,4)
        nodes_length.append('(' + nodes2[index_mini] + ': '+ str(dik) + ','+ nodes2[index_minj]+ ": "+ str(djk) + ')')
        
        # step 4
        k = nodes[-1]
        for i in dist.index : 
            if i != mini and i!=minj : 
                dist.loc[i,k] = (dist.loc[mini,i]+dist.loc[minj,i]-dist.loc[mini,minj])/2
                dist.loc[k,i] = (dist.loc[mini,i]+dist.loc[minj,i]-dist.loc[mini,minj])/2
        dist.loc[k,k]=0
        dist = dist.drop(index = [mini,minj], columns = [mini,minj])
        
    # Termination step
    i = nodes[0]
    j = nodes[1]
    m = nodes[2]
    dvi = (dist.loc[i,j]+dist.loc[i,m]-dist.loc[j,m])/2
    dvj = (dist.loc[i,j]+dist.loc[j,m]-dist.loc[i,m])/2
    dvm = (dist.loc[i,m]+dist.loc[j,m]-dist.loc[i,j])/2
    tree = "("+nodes_length[0]+":"+str(dvi)+", "+nodes_length[1]+":"+str(dvj)+", "+nodes_length[2]+":"+str(dvm)+")"
    return (tree)

########## MAIN ##########
filename = sys.argv[1]
print(NJ(filename))
