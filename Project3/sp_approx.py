import numpy as np
import sys
import random
import Bio
from Bio import SeqIO

def star_align(seqs,subst_m, gap_cost):
    #load in data
    subst_mat = read_subst_mtrx(subst_m)
    seqs_dic = read_FASTA(seqs)
    S = dic_to_list(seqs_dic)

    S1,S = find_center_seq(S, gap_cost, subst_mat)
    #print("the center seq is:")
    #print(map_back(S1,seqs_dic))
    #initialize M with the center seq:
    M = split_arr(S1,1)
    for i in range(0,len(S)):
        n = len(S1)
        m =  len(S[i])
        T = np.full([n+1,m+1], None)
        #unfortunately calculating the cost again (this has already been calculated in find_the_center_seq())
        global_linear_cost(S1,S[i],gap_cost,subst_mat,T,n,m)
        #finding the optimal alignment between the center seq and another seq
        A = IterBackTrack(S1,S[i],gap_cost,subst_mat,T)
        #print("the current pairwise to be added:")
        #print(A)
        #transposing the result, so that I work with columns rather than rows
        A = np.transpose(A)
        M = extend(M,A)    
    return M, seqs_dic
    
#helper function to translate a dictionary into a list of lists
def dic_to_list(d):
    l = []
    for seq in d.values():
        l.append(list(seq))
    return l

#helper function to translate non-nucleotides to nulceotides (ACGT)
def replace_discrepancy(seq):
    nucl = set('ACGT')
    map_nucl = {}
    seq = list(seq)
    for i in range(len(seq)):
        if seq[i] not in nucl:
            if seq[i] not in map_nucl:
                map_nucl[seq[i]] = random.sample(nucl,1)
            seq[i]=''.join(map_nucl[seq[i]])
    return ''.join(seq)

#helper function which splits an array to smaller arrays of a given size
def split_arr(arr, size):
 arrs = []
 while len(arr) > size:
     pice = arr[:size]
     arrs.append(pice)
     arr   = arr[size:]
 arrs.append(arr)
 return arrs

#this returns S1 and the rest of the seqs in a matrix
def find_center_seq(S, gap_cost, subst_mat):
    #this holds sum of scores for each seq, between itself and the other ones
    scores = np.zeros([len(S)])
    for i in range(0,len(S)):
        for j in range(i+1,len(S)):
            n = len(S[i])
            m =  len(S[j])
            L = np.full([n+1,m+1], None)
            cost = global_linear_cost(S[i],S[j],gap_cost,subst_mat,L,n,m)
            scores[i] += cost
            scores[j] += cost
    min_seq_idx = smallest(scores)
    S1 = S[min_seq_idx]
    S = [i for i in S if not np.array_equal(i, S1)]
    return S1,S

#helper function, which finds the index which holds the smallest value of that array
def smallest(arr): 
    min = arr[0] 
    min_idx = 0
    for i in range(1, len(arr)): 
        if arr[i] < min: 
            min = arr[i] 
            min_idx = i
    return min_idx

def extend(M,A):
    i=j=0
    while(i<len(M) and j<len(A)):
        if M[i][0] == '-' and A[j][0] == '-':
            M[i].append(A[j][1])
            i = i+1
            j = j+1
        elif M[i][0] == '-' and A[j][0] != '-':
            M[i].append('-')
            i = i+1
        elif M[i][0] != '-' and A[j][0] == '-':
            new_column = np.full([len(M[i])], '-', dtype=object)
            new_column = np.append(new_column,A[j][1])
            M.insert(i,new_column.tolist())
            j = j+1
            #because I actually shifted the order while inserting -> I also need to i = i+1
            i = i+1
        elif M[i][0] != '-' and A[j][0] != '-':
            M[i].append(A[j][1])
            i = i+1
            j = j+1
    #take care of the case, when you run out of one of the alignments -> add gaps
    if i<len(M):
        while i<len(M):
            M[i].append('-')
            i = i+1
    if j<len(A):
        while j<len(A):
            new_column = np.full([len(M[0])-1], '-', dtype=object)
            new_column = np.append(new_column, A[j][1])
            M.append(new_column.tolist())
            j = j+1            
    return M


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

def IterBackTrack(seq1, seq2, gap_cost, subst_mat, T):
    dict_subst = {"A":0, "C": 1, "G":2, "T":3}
    i = len(seq1)
    j = len(seq2)
    align1=align2=""

    while (i>=0 or j>=0):        
        if i>0 and j>0 and T[i,j]==T[i-1,j-1]+ subst_mat[dict_subst[seq1[i-1]], dict_subst[seq2[j-1]]] :
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i=i-1
            j=j-1
        elif i>0 and j>=0 and T[i,j]==T[i-1,j]+ gap_cost :
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i=i-1
        elif i>=0 and j>0 and T[i,j]==T[i,j-1]+gap_cost :
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j=j-1
        else: break
    return [list(align1), list(align2)]    

def read_FASTA(filename):
    records_dict = {}
    for seq_record in SeqIO.parse(filename, "fasta"):
        records_dict[seq_record.id] = replace_discrepancy(str(seq_record.seq))  
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

#helper function to extract filename without the file extension
def find_name(seqs):
    start = seqs.rfind('\\')+1
    end = seqs.find('.fasta')
    return seqs[start:end]

#helper function to find the sequence name
def map_back(seq,dic):
    seq = ''.join(seq)
    seq = seq.replace('-',"")
    inv_dic = {v: k for k, v in dic.items()}
    return inv_dic[seq]

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
    a = int(sys.argv[4])
else:
    seqs = sys.argv[1]
    subst_m = sys.argv[2]
    a = int(sys.argv[3])

alignment,seqs_dic = star_align(seqs,subst_m,5)
alignment = np.transpose(alignment)

if write_fasta:
    name = find_name(seqs)
    f = open(name+'_aligned.fasta','w')

for i in range(0,len(alignment)):
    seq_name = ">"+map_back(alignment[i],seqs_dic)
    seq_alignment = ''.join(alignment[i])
    print(seq_name)
    print(seq_alignment)
    if write_fasta:
        f.write(seq_name+" \n")
        f.write(seq_alignment+"\n")

if write_fasta:
    f.close()
