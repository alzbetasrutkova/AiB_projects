from ete3 import Tree
import sys

def get_node_names(tree):
    node_names = list()
    for node in tree:  
        node_names.append(node.name)
    return node_names


### main ###
filename1 = sys.argv[1]
filename2 = sys.argv[2]

tree1 = Tree(filename1)
tree2 = Tree(filename2)

nodes_tree1 = get_node_names(tree1)
root = nodes_tree1[0]

#print(tree1)
tree1.set_outgroup(root)
#print(tree1)
tree1.set_outgroup(root)