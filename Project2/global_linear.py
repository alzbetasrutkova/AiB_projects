import numpy as np
import sys
import math
from Bio import SeqIO

def read_FASTA(filename):
    records_dict = {}
    for seq_record in SeqIO.parse(filename, "fasta"):
        records_dict[seq_record.id] = seq_record.seq        
    return records_dict

def read_subst_mtrx(filename):
    subst_mat = np.zeros((4,4))
    f = open(filename,'r')
    f.readline()
    for i in range(0,4):
        line = f.readline()
        nums_in_line = line.split()
        for j in range(1,5):
            subst_mat[i,j-1] = nums_in_line[j]
    f.close()
    return subst_mat

def global_linear_cost(seq1, seq2, gap_cost, subst_mat, L, i, j):
    if L[i, j] == None:        
        dict_subst = {"A":0, "C": 1, "G":2, "T":3}
        v1,v2,v3,v4 = sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize
        if i>0 and j>0 : 
            v1 = global_linear_cost(seq1, seq2, gap_cost, subst_mat, L, i-1, j-1)+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]
        if i>0 and j>=0 : 
            v2 = global_linear_cost(seq1, seq2, gap_cost, subst_mat, L, i-1, j)+ gap_cost
        if i>=0 and j>0 : 
            v3 = global_linear_cost(seq1, seq2, gap_cost, subst_mat, L, i, j-1)+ gap_cost
        if i==0 and j==0 : 
            v4 = 0
        L[i,j] = min(v1,v2,v3,v4)
    return L[i,j]

def RecurBackTrack(seq1, seq2, gap_cost, subst_mat, T, i, j, align1="", align2=""):
    dict_subst = {"A":0, "C": 1, "G":2, "T":3}
        
    if i>0 and j>0 and T[i,j]==T[i-1,j-1]+ subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]] :
        align1 = seq1[i-1] + align1
        align2 = seq2[j-1] + align2
        RecurBackTrack(seq1, seq2, gap_cost, subst_mat, T,i-1, j-1, align1, align2)
    elif i>0 and j>=0 and T[i,j]==T[i-1,j]+ gap_cost :
        align1 = seq1[i-1] + align1
        align2 = "_" + align2
        RecurBackTrack(seq1, seq2, gap_cost, subst_mat, T, i-1, j, align1, align2)
    elif i>=0 and j>0 and T[i,j]==T[i,j-1]+gap_cost :
        align1 = "_" + align1
        align2 = seq2[j-1] + align2
        RecurBackTrack(seq1, seq2, gap_cost, subst_mat, T, i, j-1, align1, align2)
    else:
        f = open('alignment_with_linear_gap_cost.fasta','w')
        f.write(">Seq1 \n")
        f.write(align1)
        f.write("\n \n>Seq2 \n")
        f.write(align2)
        f.close()
        print(">Seq1")
        print (align1)
        print(">Seq2")
        print(align2)
        return None

def global_linear(seq1_filename,seq2_filename,subst_matrx_filename,gap_cost,align):
    if len(seq1_filename)>5 and seq1_filename[-5:]=="fasta" : 
        seq1 = list(read_FASTA(seq1_filename).values())[0].upper()
    else : 
        seq1 = seq1_filename.upper()
    if len(seq2_filename)>5 and seq2_filename[-5:]=="fasta":
        seq2 = list(read_FASTA(seq2_filename).values())[0].upper()
    else : 
        seq2 = seq2_filename.upper()
    subst_mat = read_subst_mtrx(subst_matrx_filename)
    
    n = len(seq1)
    m = len(seq2)
    T = np.full([n+1,m+1], None)

    C = global_linear_cost(seq1, seq2, gap_cost, subst_mat, T, n, m)
    print()
    print(C)
    print()
    
    if align:
        RecurBackTrack(seq1, seq2, gap_cost, subst_mat, T, n, m)

### main ###
align = False

if len(sys.argv) == 6:
    align = True

seq1 = ""
seq2 = ""
subst_m = ""
a = 0

if align:
    seq1 = sys.argv[2]
    seq2 = sys.argv[3]
    subst_m = sys.argv[4]
    a = int(sys.argv[5])
else:
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    subst_m = sys.argv[3]
    a = int(sys.argv[4])

global_linear(seq1, seq2,subst_m,a,align)