import sys
from ete3 import Tree

def get_node_list(t):
    node_list = []
    for node in t : 
        node_list.append(node.name)
    return node_list

def step1(t1,t2):
    t1_nodes = get_node_list(t1)
    root = t1_nodes[0]
    t1.set_outgroup(root)
    t2.set_outgroup(root)
    return t1,t2

def step2(t1,t2):
    """Assume that t1 and t2 are rooted on the same node"""
    t1_depth_first_names = []
    for node in t1.traverse("preorder"):
        if node.name != "":
            t1_depth_first_names.append(node.name)
    return (t1_depth_first_names)

def RF_dist(t1,t2):
    step1(t1,t2)
    order_t1 = step2(t1,t2) #step2
    #skip step3 because I do not use it
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
    l1 = sorted(list_intervals)
    l2 = sorted(list_intervals2)
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
    return(2*len(l1)-2*share)

### Main ###

filename1 = sys.argv[1]
filename2 = sys.argv[2]


t1 = Tree(filename1)
t2 = Tree(filename2)
print(RF_dist(t1,t2))