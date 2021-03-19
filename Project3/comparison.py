import os
import numpy as np
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

list_seq = []
for i in range(20):
    seq = "testseqs_"+str((i+1)*10)+"_3"
    list_seq.append(seq)


for i in range(20):
    file = list_seq[i]+".fasta"
    order = "python3 sp_approx.py -f " + file +" mat.txt " + str(5)
    os.system(order)
    
for i in range(20):
    file = list_seq[i]
    alignment =  list_seq[i] + "_aligned.fasta"
    order = "python3 msa_sp_score_3k.py "+ alignment
    os.system(order)
    
scores_approx = [70,135,246,328,416,491,558,648,701,727,820,910,991,1022,
                 1112,1168,1294,1243, 1458, 1494]

# Exact 

subst_mat = np.array([
        [0,5,2,5],
        [5,0,5,2],
        [2,5,0,5],
        [5,2,5,0]
    ], dtype=np.float64)

def read_FASTA(filename):
    records_dict = {}
    for seq_record in SeqIO.parse(filename, "fasta"):
        records_dict[seq_record.id] = seq_record.seq        
    return records_dict

def dic_to_list(d):
    l = []
    for seq in d.values():
        l.append(list(seq))
    return l

def C(i,j,k,a):
    if T[i,j,k] == None:        
        dict_subst = {"A":0, "C": 1, "G":2, "T":3}
        v1,v2,v3,v4,v5,v6,v7,v8 = sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize,sys.maxsize
        if i>0 and j>0 and k>0:
            v1 = C(i-1,j-1,k-1,a)+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]+subst_mat[dict_subst[seq2[j-1]], dict_subst[seq3[k-1]]]+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq3[k-1]]]
        if i>0 and j>0 and k>=0: 
            v2 = C(i-1,j-1,k,a)+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]]+(2*a)
        if i>0 and j>=0 and k>0: 
            v3 = C(i-1,j,k-1,a)+subst_mat[dict_subst[seq1[i-1]], dict_subst[seq3[k-1]]]+(2*a)
        if i>=0 and j>0 and k>0: 
            v4 = C(i,j-1,k-1,a)+subst_mat[dict_subst[seq2[j-1]], dict_subst[seq3[k-1]]]+(2*a)
        if i>0 and j>=0 and k>=0:
            v5 = C(i-1,j,k,a)+(2*a)
        if i>=0 and j>0 and k>=0:
            v6 = C(i,j-1,k,a)+(2*a)
        if i>=0 and j>=0 and k>0:
            v7 = C(i,j,k-1,a)+(2*a)
        if i==0 and j==0 and k==0:
            v8 = 0                                                                       
        T[i,j,k] = min(v1,v2,v3,v4,v5,v6,v7,v8)
    return T[i,j,k]

scores_exact = []
for i in range(20): 
    print(i)
    file = list_seq[i]+".fasta"
    seqs_dic = read_FASTA(file)
    S = dic_to_list(seqs_dic)
    seq1, seq2, seq3 = S[0], S[1], S[2]
    m = len(seq1)
    n = len(seq2)
    o = len(seq3)
    T = np.full([m+1,n+1,o+1], None)
    scores_exact.append (C(m,n,o,5))
print(scores_exact)
    
#scores_exact = [70.0, 135.0, 231.0, 318.0, 385.0, 440.0, 516.0, 589.0, 628.0, 687.0, 754.0, 810.0, 895.0, 957.0, 1023.0, 1080.0, 1186.0, 1158.0, 1323.0, 1379.0]
   
### Plot ###
xs = []
for i in range(20):
    xs.append((i+1)*10)
scores_exact = [70.0, 135.0, 231.0, 318.0, 385.0, 440.0, 516.0, 589.0, 628.0, 
                687.0, 754.0, 810.0, 895.0, 957.0, 1023.0, 1080.0, 1186.0, 
                1158.0, 1323.0, 1379.0]
scores_approx = [70,135,246,328,416,491,558,648,701,727,820,910,991,1022,
                 1112,1168,1294,1243, 1458, 1494]
ys = []
for i in range(20):
    ys.append((scores_approx[i]/scores_exact[i])/(4/3))
    
plt.plot(xs,ys)
