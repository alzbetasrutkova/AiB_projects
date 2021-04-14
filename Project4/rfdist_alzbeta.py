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


### main ###
filename1 = "C:\\Users\\alzbe\\Documents\\AU_Bioinfo_Masters\\Spring_2021\\AiB\\Projects\\Project_04\\testdata\\Testdata\\tree1.new"
filename2 = "C:\\Users\\alzbe\\Documents\\AU_Bioinfo_Masters\\Spring_2021\\AiB\\Projects\\Project_04\\testdata\\Testdata\\tree2.new"

'''
filename1 = sys.argv[1]
filename2 = sys.argv[2]
'''

t1 = Tree(filename1)
t2 = Tree(filename2)

nodes_t1 = get_node_names(t1)
root = nodes_t1[0]

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
print(df_t1)

###Step 3
num_t1 = [None] * len(df_t1) 
for i in range(0,len(df_t1)):
    num_t1[i] = i
print(num_t1)

###Step 4
map_t = map_arrs(df_t1,num_t1)

print(map_t)

#print(t1)
#print(t1.get_common_ancestor("seq6", "seq1"))

intervals_t1 = get_df_intervals(t1, map_t)
intervals_t2 = get_df_intervals(t2, map_t)

print(intervals_t1)
print(intervals_t2)

###Step 5
to_be_sorted = intervals_t1 + intervals_t2
print(to_be_sorted)
'''

###Step 4
intervals_t1 = step4(t1)
intervals_t2 = step4(t2)

print(intervals_t1)
print(intervals_t2)

'''


