import sys
from ete3 import Tree

#We are using the ete3 package to work with trees
#The steps in this implementation correspond to the algorithm steps as presented in slides "Tree_comparison.pdf" (23-27)

#function which returns a list of all the node names in a given tree
def get_node_list(t):
    node_list = []
    for node in t : 
        node_list.append(node.name)
    return node_list

#rooting the trees
#the root is chosen as the first node in tree 1
#we have also tried rooting with a fixed species, which didn't help 
def step1(t1,t2):
    t1_nodes = get_node_list(t1)
    root = t1_nodes[0]
    #root = "4__gi|8134331|sp|Q9Y2Q0.1|AT8A1_HUMAN"
    t1.set_outgroup(root)
    t2.set_outgroup(root)
    return t1,t2

#finding the DFS numbergin for the first tree
def step2(t1,t2):
    t1_depth_first_names = []
    for node in t1.traverse("preorder"):
        if node.name != "":
            t1_depth_first_names.append(node.name)
    return (t1_depth_first_names)

#a wrapper function
def RF_dist(t1,t2):
    step1(t1,t2)
    order_t1 = step2(t1,t2) 
    #step 3 is done in steps 4
    #step 4.1
    t1_nodes_children = []
    for node in t1.traverse("preorder"):
        if node.name == "":
            t1_nodes_children.append(node.get_leaf_names())
    list_intervals = []
    for l in t1_nodes_children:
        m = order_t1.index(l[0])
        M = order_t1.index(l[0])
        for i in range(len(l)):
            if order_t1.index(l[i])<m:
                m= order_t1.index(l[i])
            if order_t1.index(l[i])>M : 
                M= order_t1.index(l[i])
        list_intervals.append("["+str(m)+","+str(M)+"]")
    
    # Step 4.2 
    t2_nodes_children = []
    for node in t2.traverse("preorder"):
        if node.name == "":
            t2_nodes_children.append(node.get_leaf_names())
    list_intervals2 = []
    for l in t2_nodes_children:
        order_t1 = step2(t1,t2)
        m = order_t1.index(l[0])
        M = order_t1.index(l[0])
        for i in range(len(l)):
            if order_t1.index(l[i])<m:
                m= order_t1.index(l[i])
            if order_t1.index(l[i])>M : 
                M= order_t1.index(l[i])
        if M-m+1 == len(l):
            list_intervals2.append("["+str(m)+","+str(M)+"]")

    # Step 5 
    #print("Lenght of the first interval:", len(list_intervals))
    #print(list_intervals)
    #print("Lenght of the second interval:", len(list_intervals2))
    #print(list_intervals2)

    l1 = sorted(list_intervals)
    l2 = sorted(list_intervals2)
    '''
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
    '''
    share = len(list(set(list_intervals) & set(list_intervals2)))
    #print("Num of shared intervals:",share)
    return(2*len(l1)-2*share)

### Main ###

filename1 = sys.argv[1]
filename2 = sys.argv[2]


t1 = Tree(filename1)
t2 = Tree(filename2)

print(RF_dist(t1,t2))