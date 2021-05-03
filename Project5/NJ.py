# -*- coding: utf-8 -*-
"""
Created on Sun May  2 10:04:59 2021

@author: lenab
"""

import numpy as np
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
    S = dist.shape[1]
    N = np.zeros((S,S))
    rs = []
  
    for i in range(S):
        ri = 0
        for j in range(S):
            ri +=dist[i,j]
        rs.append(ri/(S-2))
          
    minN = sys.maxsize
    index_mini = 0
    index_minj = 1
   
    for i in range(S):
        for j in range(S):
            nij = dist[i,j]-rs[i]-rs[j]
            N[i,j]=nij
            if i != j:
                if nij<minN : 
                    minN = nij
                    index_mini = i
                    index_minj = j
    rimin = rs[index_mini]
    rjmin = rs[index_minj]
    return(index_mini, index_minj,rimin,rjmin)


def NJ(dist_mat):
    # create the distance dataframe
    species_dict, dist = read_mtrx(dist_mat)
    species = list(species_dict.keys())
    # initializing variables
    nodes_length = species.copy()
    nodes = species.copy()
    
    while len(nodes)>3 : 
        
        #step 1 
        index_mini, index_minj, rimin, rjmin = step1(dist)
        
        # step 2 
        #nodes = list(dist.index).copy()
        nodes_copy = nodes.copy()
        mini = nodes[index_mini]
        minj = nodes[index_minj]
        nodes.remove(mini)
        nodes.remove(minj)
        nodes.append('('+mini+','+minj+')')
        
        #step 3 
        nodes2 = nodes_length.copy()
        nodes_length.remove(nodes2[index_mini])
        nodes_length.remove(nodes2[index_minj])
        dik = round((dist[index_mini,index_minj] + rimin - rjmin)/2,4)
        djk = round((dist[index_mini,index_minj] + rjmin - rimin)/2,4)
        nodes_length.append('(' + nodes2[index_mini] + ': '+ str(dik) + ','+ nodes2[index_minj]+ ": "+ str(djk) + ')')
      
        # step 4
        newrow = []
        for i in range(len(nodes_copy)) :
            dik2 = (dist[index_mini,i]+dist[index_minj,i]-dist[index_mini,index_minj])/2
            newrow.append(dik2)
        dist = np.vstack([dist,newrow])
        newrow.append(0)
        newrow = np.array((newrow))
        dist = np.c_[dist,newrow]
        dist = np.delete(dist, np.s_[index_mini, index_minj], axis=0)
        dist = np.delete(dist, np.s_[index_mini, index_minj], axis=1)
    
    # Termination step
    dvi = (dist[0,1]+dist[0,2]-dist[1,2])/2
    dvj = (dist[0,1]+dist[1,2]-dist[0,2])/2
    dvm = (dist[0,2]+dist[1,2]-dist[0,1])/2
    tree = "("+nodes_length[0]+":"+str(dvi)+", "+nodes_length[1]+":"+str(dvj)+", "+nodes_length[2]+":"+str(dvm)+")"
    return (tree)

########## MAIN ##########

write_fasta = False

if len(sys.argv) == 2:
    write_fasta = True
    filename = sys.argv[2]
else : 
    filename = sys.argv[1]
    
tree = NJ(filename)

if write_fasta==True:
    name = filename[:-4]
    f = open(name + '_NJtree.newick','w')
    f.write(tree)
    f.close()
    

