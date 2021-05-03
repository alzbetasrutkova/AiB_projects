from ete3 import Tree
import numpy as np
import sys

def get_node_names(t):
    node_names = list()
    for node in t:  
        node_names.append(node.name)
    return node_names

def dfs_numbering(t):
    depth_first_names = []
    for node in t.traverse("preorder"):
        if node.name != "":
            depth_first_names.append(node.name)
    return(depth_first_names)

def map_arrs(arr1,arr2):
    dic = {}
    i = 0
    for k in arr1:
        dic[k] = arr2[i]
        i = i+1
    return dic

def is_it_interval(min,max,size):
    if (max-min+1 == size):
        return True
    else:
        return False

#this finds the intervals, but runs in O(n^2) not O(n) as required
def get_df_intervals(t, map_t):
    print(map_t)
    inters = list()
    for node in t.traverse("preorder"):
        if node.name == "":
            min = sys.maxsize
            max = 0
            s = 0
            for n in node.traverse():
                if n.name != "":
                    #print(map_t[n.name])
                    if (map_t[n.name] < min):
                        min = map_t[n.name]
                    if (map_t[n.name] > max):
                        max = map_t[n.name]
                    s = s+1
            if is_it_interval(min,max,s):
                inter = [min, max]
                inters.append(inter)            
    return inters    

def find_shared(arr):
    sorted_arr = sorted(arr)
    #print(sorted_arr)
    i = 0
    shared = 0
    while i <=len(arr):
        if i + 1 < len(arr):
            if sorted_arr[i] == sorted_arr[i+1]:
                shared = shared + 1
                i = i + 2
            else:
                i = i + 1
        else:
            i = i + 1
    return shared



### main ###
'''
filename1 = "C:\\Users\\alzbe\\Documents\\AU_Bioinfo_Masters\\Spring_2021\\AiB\\Projects\\Project_04\\testdata\\Testdata\\tree1.new"
filename2 = "C:\\Users\\alzbe\\Documents\\AU_Bioinfo_Masters\\Spring_2021\\AiB\\Projects\\Project_04\\testdata\\Testdata\\tree2.new"

'''
filename1 = sys.argv[1]
filename2 = sys.argv[2]


t1 = Tree(filename1)
t2 = Tree(filename2)

nodes_t1 = get_node_names(t1)
root = nodes_t1[0]
#root = "4__gi|8134331|sp|Q9Y2Q0.1|AT8A1_HUMAN"

###Step 1
#print("Before rooting:")
#print(t1)
#print(t2)

#print("After rooting:")
t1.set_outgroup(root)
#print(t1)
t2.set_outgroup(root)
#print(t2)

###Step 2
df_t1 = dfs_numbering(t1)
#print(df_t1)

###Step 3
num_t1 = [None] * len(df_t1) 
for i in range(0,len(df_t1)):
    num_t1[i] = i
#print(num_t1)

###Step 4
map_t = map_arrs(df_t1,num_t1)

#print(map_t)

#print(t1)
#print(t1.get_common_ancestor("seq6", "seq1"))

intervals_t1 = get_df_intervals(t1, map_t)
intervals_t2 = get_df_intervals(t2, map_t)

print("The length of the first intervals:", len(intervals_t1))
print("The length of the second intervals:", len(intervals_t2))


###Step 5
pool = intervals_t1 + intervals_t2
#print(pool)

num_shared = find_shared(pool)
#print(len(pool))
print("Num of shared intervals", num_shared)

print("The RF distance is: ", (len(pool) - (num_shared*2))*2)
'''
l1 = sorted(intervals_t1)
l2 = sorted(intervals_t2)
i=0
j=0
share = 0
while i<len(l1) and j<len(l2):
    if l1[i]==l2[j]:
        share+=1
        i+=1
        j+=1
    else : 
        i+=1

print(2*len(l1)-2*share)
'''


