import math 
import sys
import Bio
from Bio import SeqIO
import numpy as np

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

def cost_affine(seq1,seq2,subst_mat,a,b):
    m = len(seq1)
    n = len(seq2)
    
    S = np.full([m+1,n+1], None)
    I = np.full([m+1,n+1], None)
    D = np.full([m+1,n+1], None)
    
    dict_subst = {"A":0, "C": 1, "G":2, "T":3}
    
    for i in range(0, m+1):
        for j in range(0, n+1):
        # Compute D[i,j]
            v1,v2 = sys.maxsize,sys.maxsize
            if i>0 and j>=0 : 
                v1 = S[i-1,j]+(a+b)
            if i>1 and j>=0:
                v2 = D[i-1,j]+a
            D[i,j] = min(v1,v2)
                    
        # Compute I[i,j]
            v1,v2 = sys.maxsize,sys.maxsize
            if i>=0 and j>0:
                v1 = S[i,j-1]+(a+b)
            if i>=0 and j>1:
                v2 = I[i,j-1]+a
            I[i,j] = min(v1,v2)
                 
        # Compute S[i,j] 
            v1,v2,v3,v4 = sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize
            if i==0 and j==0 :
                v1 = 0                
            if i>0 and j>0:
                v2 = S[i-1,j-1] + subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]
            if i>0 and j>=0:
                v3 = D[i,j]
            if i>=0 and j>0:
                v4 = I[i,j]
            S[i,j] = min(v1,v2,v3,v4)        
    return S

def backtrack_affine(seq1,seq2,S,subst_mat,a,b):
    i = len(seq1)
    j = len(seq2)
    align1 = ""
    align2 = ""
    dict_subst = {"A":0, "C": 1, "G":2, "T":3}
    while (i>0 or j>0):
        if (i>0 and j>0) and (S[i,j] == S[i-1,j-1] + subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]):
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i = i-1
            j = j-1
        else:
            k = 1
            while True:
                if i>=k and S[i,j] == S[i-k,j] + (a+b*k):
                    l = i
                    while(l>=i-k+1):
                        align1 = seq1[l-1] + align1
                        align2 = "-" + align2
                        l = l-1
                    i = i-k
                    break
                elif j>=k and S[i,j] == S[i,j-k] + (a+b*k):
                    l = j
                    while(l>=i-k+1):
                        align1 = "-" + align1
                        align2 = seq2[l-1] + align2
                        l = l-1
                    j = j-k
                    break
                else:
                    k = k+1
    return align1, align2

def FASTA_out(align1,align2):
    print(">seq1")
    print(align1)
    print(">seq2")
    print(align2) 

def global_affine(seq1_filename,seq2_filename,subst_matrx_filename,a,b,align):
    if len(seq1_filename)>5 and seq1_filename[-5:]=="fasta" : 
        seq1 = list(read_FASTA(seq1_filename).values())[0].upper()
    else: 
        seq1 = seq1_filename.upper()
    if len(seq2_filename)>5 and seq2_filename[-5:]=="fasta":
        seq2 = list(read_FASTA(seq2_filename).values())[0].upper()
    else: 
        seq2 = seq2_filename.upper()
    
    subst_mat = read_subst_mtrx(subst_matrx_filename)

    S = cost_affine(seq1,seq2,subst_mat,a,b)
    print()
    print(S[len(seq1),len(seq2)])
    print()
    
    if align:
        align1,align2 = backtrack_affine(seq1,seq2,S,subst_mat,a,b)
        FASTA_out(align1,align2)

### main ### 
align = False

if len(sys.argv) == 7:
    align = True

seq1 = ""
seq2 = ""
subst_m = ""
a = 0
b = 0

if align:
    seq1 = sys.argv[2]
    seq2 = sys.argv[3]
    subst_m = sys.argv[4]
    a = int(sys.argv[5])
    b = int(sys.argv[6])
else:
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    subst_m = sys.argv[3]
    a = int(sys.argv[4])
    b = int(sys.argv[5])

global_affine(seq1, seq2,subst_m,a,b,align)



