import sys
import numpy as np
from ete3 import Tree


#an improved version of the matrix reader for project_3
def read_subst_mtrx(filename):
    f = open(filename,'r')
    num_species = int(f.readline())

    dict_mat = {}
    mat = np.zeros((num_species,num_species))
    
    for i in range(0,num_species):
        line = f.readline()
        nums_in_line = line.split()
        dict_mat[nums_in_line[0]] = i
        for j in range(1,num_species):
            mat[i,j-1] = nums_in_line[j]
    f.close()
    return dict_mat, mat 


def tree_contruct_nj(dict_subst, mat ):
    start = True

    while(dict_mat.size > 3):
        if start:
            mat = 
            start = False
            pass
        else:


### main ###
matr_filename = "C:\Users\awila\Downloads\example_slide4.phy"

print(read_subst_mtrx(matr_filename))