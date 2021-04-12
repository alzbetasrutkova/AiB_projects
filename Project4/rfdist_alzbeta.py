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

def get_df_intervals(t, df_num):
    inters = list()
    for node in t.traverse("preorder"):
        if node.name == "":
            min = 0, max = 0, size = 0
            for n in node.traverse():
                if 
            
    return(inters)


### main ###
filename1 = sys.argv[1]
filename2 = sys.argv[2]

tree1 = Tree(filename1)
tree2 = Tree(filename2)

nodes_tree1 = get_node_names(tree1)
root = nodes_tree1[0]

###Root the trees
#print(tree1)
tree1.set_outgroup(root)
#print(tree1)
tree1.set_outgroup(root)

###Get the DFS numbering
df_tree1 = dfs_numbering(tree1)
#print(df_tree1)
df_tree2 = dfs_numbering(tree2)
#print(df_tree2)

###Get the DF-intervals







