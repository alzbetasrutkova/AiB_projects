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

def C(seq1, seq2, seq3, a, subst_mat, T, i, j,k):
    if T[i,j,k] == None:        
        dict_subst = {"A":0, "C": 1, "G":2, "T":3}
        v1,v2,v3,v4,v5,v6,v7,v8 = sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize
        if i>0 and j>0 and k>0:
            v1 = C(seq1,seq2,seq3,a,subst_mat,T,i-1,j-1,k-1)+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]+subst_mat[dict_subst[seq2[j-1]], dict_subst[seq3[k-1]]]+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq3[k-1]]]
        if i>0 and j>0 and k>=0: 
            v2 = C(seq1,seq2,seq3,a,subst_mat,T,i-1,j-1,k)+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]+(2*a)
        if i>0 and j>=0 and k>0: 
            v3 = C(seq1,seq2,seq3,a,subst_mat,T,i-1,j,k-1)+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq3[k-1]]]+(2*a)
        if i>=0 and j>0 and k>=0: 
            v4 = C(seq1,seq2,seq3,a,subst_mat,T,i,j-1,k-1)+subst_mat[dict_subst[seq2[j-1]], dict_subst[seq3[k-1]]]+(2*a)
        if i>0 and j>=0 and k>=0:
            v5 = C(seq1,seq2,seq3,a,subst_mat,T,i-1,j,k)+(2*a)
        if i>=0 and j>0 and k>=0:
            v6 = C(seq1,seq2,seq3,a,subst_mat,T,i,j-1,k)+(2*a)
        if i>=0 and j>=0 and k>0:
            v7 = C(seq1,seq2,seq3,a,subst_mat,T,i,j,k-1)+(2*a)
        if i==0 and j==0 and k==0:
            v8 = 0                                                                       
        T[i,j,k] = min(v1,v2,v3,v4,v5,v6,v7,v8)
    return T[i,j,k]

def RecurBackTrack_write(seq1, seq2, seq3, a, subst_mat, T, i,j,k, align1="", align2="", align3=""):
    dict_subst = {"A":0, "C": 1, "G":2, "T":3}
        
    if i>0 and j>0 and k>0 and T[i,j,k]==T[i-1,j-1, k-1]+ subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]+subst_mat[dict_subst[seq2[j-1]], dict_subst[seq3[k-1]]]+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq3[k-1]]]:
        align1 = seq1[i-1] + align1
        align2 = seq2[j-1] + align2
        align3 = seq3[k-1] + align3
        RecurBackTrack(seq1, seq2, seq3,  a, subst_mat, T,i-1, j-1, k-1, align1, align2, align3)
    elif i>0 and j>0 and k>=0 and T[i,j,k]==T[i-1,j-1,k]+ subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]+(2*a):
        align1 = seq1[i-1] + align1
        align2 = seq2[j-1] + align2
        align3 = "_" + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i-1, j-1, k, align1, align2, align3)
    elif i>0 and j>=0 and k>0 and T[i,j,k]==T[i-1,j, k-1]+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq3[k-1]]]+(2*a) :
        align1 = seq1[i-1]+ align1
        align2 = "_" + align2
        align3 = seq3[k-1] + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i-1, j, k-1, align1, align2, align3)
    elif i>=0 and j>0 and k>0 and T[i,j,k]==T[i, j-1, k-1] + subst_mat[dict_subst[seq2[j-1]], dict_subst[seq3[k-1]]]+(2*a) :
        align1 = "_" + align1
        align2 = seq2[j-1] + align2
        align3 = seq3[k-1] + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i, j-1, k-1, align1, align2, align3)
    elif i>0 and j>=0 and k>=0 and T[i,j,k] == T[i-1, j, k] + (2*a) : 
        align1 = seq1[i-1] + align1
        align2 = "_" + align2
        align3 = "_" + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i-1, j, k, align1, align2, align3)
    elif i>=0 and j>0 and k>=0 and T[i,j,k] == T[i,j-1,k] + (2*a) : 
        align1 = "_" + align1
        align2 = seq2[j-1] + align2
        align3 = "_" + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i, j-1, k, align1, align2, align3)
    elif i>=0 and j>=0 and k>0 and T[i,j,k]== T[i,j,k-1] + (2*a) : 
        align1 = "_" + align1
        align2 = "_" + align2
        align3 = seq3[k-1]+align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i, j, k-1, align1, align2, align3)
    else : 
        f = open('alignment_with_linear_gap_cost.fasta','w')
        f.write(">Seq1 \n")
        f.write(align1)
        f.write("\n \n>Seq2 \n")
        f.write(align2)
        f.write("\n \n>Seq3 \n")
        f.write(align3)
        f.close()
        print(">Seq1")
        print (align1)
        print(">Seq2")
        print(align2)
        print(">Seq3")
        print(align3)
        return None

def RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i,j,k, align1="", align2="", align3=""):
    dict_subst = {"A":0, "C": 1, "G":2, "T":3}
    print(type(dict_subst[seq1[i-1]]))
    print(j)
    print(k)
        
    if i>0 and j>0 and k>0 and T[i,j,k]==T[i-1,j-1, k-1]+ subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]+subst_mat[dict_subst[seq2[j-1]], dict_subst[seq3[k-1]]]+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq3[k-1]]]:
        align1 = seq1[i-1] + align1
        align2 = seq2[j-1] + align2
        align3 = seq3[k-1] + align3
        RecurBackTrack(seq1, seq2, seq3,  a, subst_mat, T,i-1, j-1, k-1, align1, align2, align3)
    elif i>0 and j>0 and k>=0 and T[i,j,k]==T[i-1,j-1,k]+ subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]+(2*a):
        align1 = seq1[i-1] + align1
        align2 = seq2[j-1] + align2
        align3 = "_" + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i-1, j-1, k, align1, align2, align3)
    elif i>0 and j>=0 and k>0 and T[i,j,k]==T[i-1,j, k-1]+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq3[k-1]]]+(2*a) :
        align1 = seq1[i-1]+ align1
        align2 = "_" + align2
        align3 = seq3[k-1] + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i-1, j, k-1, align1, align2, align3)
    elif i>=0 and j>0 and k>0 and T[i,j,k]==T[i, j-1, k-1] + subst_mat[dict_subst[seq2[j-1]], dict_subst[seq3[k-1]]]+(2*a) :
        align1 = "_" + align1
        align2 = seq2[j-1] + align2
        align3 = seq3[k-1] + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i, j-1, k-1, align1, align2, align3)
    elif i>0 and j>=0 and k>=0 and T[i,j,k] == T[i-1, j, k] + (2*a) : 
        align1 = seq1[i-1] + align1
        align2 = "_" + align2
        align3 = "_" + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i-1, j, k, align1, align2, align3)
    elif i>=0 and j>0 and k>=0 and T[i,j,k] == T[i,j-1,k] + (2*a) : 
        align1 = "_" + align1
        align2 = seq2[j-1] + align2
        align3 = "_" + align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i, j-1, k, align1, align2, align3)
    elif i>=0 and j>=0 and k>0 and T[i,j,k]== T[i,j,k-1] + (2*a) : 
        align1 = "_" + align1
        align2 = "_" + align2
        align3 = seq3[k-1]+align3
        RecurBackTrack(seq1, seq2, seq3, a, subst_mat, T, i, j, k-1, align1, align2, align3)
    else : 
        print(">Seq1")
        print (align1)
        print(">Seq2")
        print(align2)
        print(">Seq3")
        print(align3)
        return None

def sp_exact(seqs,subst_matrx_filename,gap_cost,write_fasta):
    seq1 = list(read_FASTA(seqs).values())[0].upper()
    seq2 = list(read_FASTA(seqs).values())[1].upper()
    seq3 = list(read_FASTA(seqs).values())[2].upper()

    m = len(seq1)
    n = len(seq2)
    o = len(seq3)

    T = np.full([m+1,n+1,o+1], None)

    subst_mat = read_subst_mtrx(subst_matrx_filename)

    cost = C(seq1, seq2, seq3, gap_cost, subst_mat, T, m, n, o)
    print(cost)

    if write_fasta:
        RecurBackTrack_write(seq1,seq2,seq3,gap_cost,subst_m,T,m,n,o)
    else:
        RecurBackTrack(seq1,seq2,seq3,gap_cost,subst_m,T,m,n,o)

### main ###
write_fasta = False

if len(sys.argv) == 5:
    write_fasta = True

seqs = ""
subst_m = ""
gap_cost = 0

if write_fasta:
    seqs = sys.argv[2]
    subst_m = sys.argv[3]
    gap_cost = int(sys.argv[4])
else:
    seqs = sys.argv[1]
    subst_m = sys.argv[2]
    gap_cost = int(sys.argv[3])

sp_exact(seqs,subst_m,gap_cost,write_fasta)
