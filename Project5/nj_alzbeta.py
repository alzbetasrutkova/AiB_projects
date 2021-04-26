import sys
import numpy as np
from ete3 import Tree


#an improved version of the read_subst_mtrx from project_3
#returns the matrix and a dictionary mapping the species name to indeces
def read_mtrx(filename):
    f = open(filename,'r')
    num_species = int(f.readline())

    dict_mat = {}
    D = np.zeros((num_species,num_species))
    
    for i in range(0,num_species):
        line = f.readline()
        nums_in_line = line.split()
        dict_mat[nums_in_line[0]] = i
        for j in range(1,num_species +1):
            D[i,j-1] = nums_in_line[j]
    f.close()
    return dict_mat, D 

#don't have to search the whole matrix, because it is symetric
#also no need to look at the diagonal?
def find_min(N):
    min_i = 0
    min_j = 0

    num_rows, num_cols = N.shape
    counter = 1
    min_val = sys.maxsize

    for i in range(1,num_rows):
        for j in range(0,counter):
            if min_val > N[i,j]:
                min_val = N[i,j]
                min_i = i
                min_j = j                
        counter = counter+1

    return min_i, min_j

'''
#NOT CURRENTLY USING 
#might not work for a matrix which is not n*n
def compute_r(D,index,row):
    r = 0
    dim_D = D.shape[1]
    if row:
        for i in range(0,dim_D):
            r = r+D[index, i]
    else:
        for i in range(0,dim_D):
            r = r+D[i, index]
    return r*(1/(dim_D-2))
'''

def find_ris(D):
    ris = []
    dim_D = D.shape[1]
    for i in range(dim_D):
        ri = 0
        for j in range(dim_D):
            ri +=D[i,j]
        ris.append(ri/(dim_D-2))
    return ris

def find_rjs(D):
    rjs = []
    dim_D = D.shape[1]
    for j in range(dim_D):
        rj = 0
        for i in range(dim_D):
            rj +=D[i,j]
        rjs.append(rj/(dim_D-2))
    return rjs

def compute_N(D):
    num_rows, num_cols = D.shape
    N = np.zeros((num_rows,num_cols))

    ris = find_ris(D)
    rjs = find_rjs(D)

    for i in range(0,num_rows):
        for j in range(0, num_cols):
            #r_i = compute_r(D, i, True)
            #r_j = compute_r(D, j, False)
            N[i,j] = D[i,j] - (ris[i] + rjs[j])
    return N

def find_key_by_val(d, val):
    for v, k in d.items():  
        if k == val:
            return v

def grow_tree(dict_species, D, min_i, min_j, t_newick):
    ris = find_ris(D)
    rjs = find_rjs(D)

    #print(dict_species)
    spec_i = find_key_by_val(dict_species, min_i)
    spec_j = find_key_by_val(dict_species, min_j)

    weight_i = 0.5*(D[min_i, min_j] + ris[min_i] - rjs[min_j])
    weight_j = D[min_i, min_j] - weight_i

    t_newick = "("+spec_i+":"+str(weight_i)+","+spec_j+":"+str(weight_j)+")"

    return spec_i, spec_j, t_newick
     
#not finished yet
def update_D(D, min_i, min_j):
    #first delete the corresponding rows and columns
    new_D = np.delete(D, min_i, min_j)
    new_D = np.delete(new_D, min_j, min_i)

    

    return new_D

def tree_construct_nj(dict_species, D, t_newick ):

    while(D.size > 3):
        #step 1
        N = compute_N(D)
        #print(N)

        #step 2
        min_i, min_j = find_min(N)
        #print(min_i, min_j)

        #step 3
        #add a new node = add new set of brackets? ( , )
        spec_i, spec_j, t_newick = grow_tree(dict_species, D, min_i, min_j, t_newick)
        print(t_newick)

        #step 4
        D = update_D(D, min_i, min_j)
        break


### main ###
matr_filename = "C:\\Users\\alzbe\\Documents\\AU_Bioinfo_Masters\\Spring_2021\\AiB\\Projects\\Project_05\\example_slide4.phy"

dict_species, D = read_mtrx(matr_filename)

t_newick = ""

tree_construct_nj(dict_species, D, t_newick)
