import sys
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
#seq_in = "hphhhhhh"
seq_in = "hhppppphhppphppphp"
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

print(pairs)
print(score_left)
print(score_right)

#print(pairs[0][0])

directions = ""

#S1
if pairs[0][0] != 0:
    i = pairs[0][0]
    while i > 0:
        directions = directions + "f"
        i = i-1      

for i in range(0, len(pairs)):
    if i < len(pairs)-1:
        between = pairs[i+1][0] - pairs[i][0] - 1   
        if between == 1:
            directions = directions + "ff"
        else:
            #there is always an odd num of elements between two matches        
            left_part = (between / 2) - 1
            j = left_part
            directions = directions + "l"
            while(j > 0):
                directions = directions + "f"
                j = j - 1
            directions = directions + "r"
            j = left_part
            while(j > 0):
                directions = directions + "f"
                j = j - 1

#creating the main fold:
between_s1_s2 = pairs[len(pairs)-1][1] - pairs[len(pairs)-1][0] - 1
j = between_s1_s2 / 2
while(j > 0):
    directions = directions + "f"
    j = j-1
directions = directions + "rr"
j = (between_s1_s2 / 2) - 1
while(j > 0):
    directions = directions + "f"
    j = j-1

#S2
for i in range(len(pairs)-1, 0, -1):
    if i > 0:
        between = pairs[i-1][1] - pairs[i][1] - 1  
        if between == 1:
            directions = directions + "ff"
        else:
            #there is always an odd num of elements between two matches        
            left_part = (between / 2) - 1
            j = left_part
            directions = directions + "l"
            while(j > 0):
                directions = directions + "f"
                j = j - 1
            directions = directions + "r"
            j = left_part
            while(j > 0):
                directions = directions + "f"
                j = j - 1
print(seq_in)
print(directions)
             











