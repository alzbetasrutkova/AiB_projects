import sys
import numpy as np
from ete3 import Tree


#an improved version of the read_subst_mtrx from project_3
#returns the matrix and a dictionary mapping the species name to indeces
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

def compute_mat():
    pass

def tree_contruct_nj(dict_subst, mat ):
    start = True

    while(dict_mat.size > 3):
        if not start:
            mat = compute_mat()
            start = False
            pass
        else:
            pass


### main ###
matr_filename = "C:\Users\awila\Downloads\example_slide4.phy"

print(read_mtrx(matr_filename))