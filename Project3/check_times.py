# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 15:22:04 2021

@author: lenab
"""
 
import os
import time

seqs = ["3seqs.fasta", "4seqs.fasta", "5seqs.fasta", "6seqs.fasta", "brca1-testseqs.fasta"]

times = []
for i in range(len(seqs)):
    file = seqs[i]
    order = "python3 sp_approx.py -f " + file + " mat.txt 5"
    start_time = time.time()
    os.system(order)
    times.append(time.time() - start_time)
print (times)
# [0.9413280487060547, 1.4405720233917236, 2.0659964084625244, 2.7901368141174316, 4.654843091964722]

order = "python3 sp_exact_3.py -f 3seqs.fasta mat.txt 5"
start_time = time.time()
os.system(order)
print(time.time()-start_time)
# 86.27345585823059