from Bio import Phylo
import sys

### main ###
filename1 = sys.argv[1]
filename2 = sys.argv[2]

tree1 = Phylo.read(filename1, "newick")
tree2 = Phylo.read(filename2, "newick")

Phylo.draw_ascii(tree1)

tree1.root_with_outgroup({"name": "seq5"})
Phylo.draw_ascii(tree1)


