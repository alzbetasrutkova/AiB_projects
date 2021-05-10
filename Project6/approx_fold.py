import sys
import numpy as np
'''
Implement the 1/4 approximation algorithm for folding a hp-string in the 2D HP Model. You program should take a hp-string as input
 and output the energy of the computed fold for this string. It should also be possible to output a representation of the computed
  fold in a format that can be used by hpview.py (or hpview3k.py if you are using Python 3).
'''

def odds_evens(seq):
    odd_even = [None] * len(seq)
    for i in range(0, len(seq)):
        if seq[i] == "h":
            if i % 2 == 0:
                #an even h
                odd_even[i] = "e"
            else:
                #an odd h
                odd_even[i] = "o"
        else:
            #a p
            odd_even[i] = "p"
    return odd_even


def find_score_and_pairs(seq_odd_even, evens_left = True):
    score = 0
    pairs = list()
    j = len(seq_odd_even) - 1
    i = 0

    if evens_left:
        while i < j:
            if seq_odd_even[i] == "e":
                while j > i:
                    if seq_odd_even[j] == "o":
                        score = score + 1
                        pairs.append([i, j])
                        j = j-1
                        break
                    j = j-1
            i = i + 1
    else:
        while i < j:
            if seq_odd_even[i] == "o":
                while j > i:
                    if seq_odd_even[j] == "e":
                        score = score + 1
                        pairs.append([i, j])
                        j = j-1
                        break
                    j = j-1
            i = i + 1
    return score, pairs

                




### main ###
#seq_in = sys.argv[1]

#assign the odd/even labels
seq_in = "hphhhhhh"
seq_odd_even = odds_evens(seq_in)
#print(seq_odd_even)

#find the scores and pairs for starting with evens/odds
score_left, pairs_left = find_score_and_pairs(seq_odd_even, True)
score_right, pairs_right = find_score_and_pairs(seq_odd_even, False)

#choose the max
pairs = list()
if score_left >= score_right:
    pairs = pairs_left
else:
    pairs = pairs_right





