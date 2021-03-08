import numpy as np
import sys

def star_align(S, gap_cost,subst_mat):
    S1,S = find_center_seq(S, gap_cost, subst_mat)
    #print(S1)
    #print(S)
    #initialize M with the center seq:
    M = split_arr(S1,1)
    print(M)    
    #print(M[0][0])
    n = len(S1)
    m =  len(S[0])
    T = np.full([n+1,m+1], None)
    #unfortunately calculating the cost again:/
    global_linear_cost(S1,S[0],gap_cost,subst_mat,T,n,m)
    #finding the optimal alignment between the center seq and another seq
    A = IterBackTrack(S1,S[0],gap_cost,subst_mat,T)
    #transposing the result, so that I work with columns rather than rows
    A = np.transpose(A)
    #print(A[1][0])
    #print(A.tranpose)
    M = extend(M,A)
    print(M)
    '''
    for i in range(0,len(S)):
        A = RecurBackTrack(S1,S[i])
        M = extend(M,A)
    return M
    '''

#helper function which splits an array arr to smaller arrays of a given size
def split_arr(arr, size):
 arrs = []
 while len(arr) > size:
     pice = arr[:size]
     arrs.append(pice)
     arr   = arr[size:]
 arrs.append(arr)
 return arrs

#this will return S1 and the rest of the seqs in a matrix
def find_center_seq(S, gap_cost, subst_mat):
    #this will hold sum of scores for each seq, between itself and the other ones
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

#helper function, which finds the index of an array which holds the smallest value of that array
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
    #take care of the case, when you run out of one of the alignments -> add gaps
    while(i<len(M) and j<len(A)):
        if M[i][0] == '_' and A[j][0] == '_':
            M[i].append(A[j][1])
            i = i+1
            j = j+1
        elif M[i][0] == '_' and A[j][0] != '_':
            M[i].append('_')
            i = i+1
        elif M[i][0] != '_' and A[j][0] == '_':
            gap_column = np.full([len(M[i])], '_')
            M.insert(i,gap_column)
            j = j+1
        elif M[i][0] != '_' and A[j][0] != '_':
            M[i].append(A[j][1])
            i = i+1
            j = j+1
    if i<len(M)
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
            align2 = "_" + align2
            i=i-1
        elif i>=0 and j>0 and T[i,j]==T[i,j-1]+gap_cost :
            align1 = "_" + align1
            align2 = seq2[j-1] + align2
            j=j-1
        else: break
    return [list(align1), list(align2)]    


### main ###

S = [['A','A','C','T','G'], 
    ['A','C','C','T'],
    ['A','C','T']]

subst_mat = np.array([
        [0,5,2,5],
        [5,0,5,2],
        [2,5,0,5],
        [5,2,5,0]
    ])

#S1, S_rest = find_center_seq(S,5,subst_mat)
#print(S)
#print(S1)
#print(S_rest)
star_align(S,5,subst_mat)
#L = np.full([n+1,m+1], None)

