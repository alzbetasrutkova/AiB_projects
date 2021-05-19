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


##################Lena's code#######################
def score(seq, fold):
    if seq[0]=="h" or seq[0]=="H":
        grid={(0,0):0}
    else : 
        grid = {}
    i=0
    j=0
    flag = "right"
    for k in range(len(fold)):
        if fold[k]=="f" : 
            if flag == "right":
                i+=1
            elif flag=="left":
                i-=1
            elif flag == "up":
                j+=1
            else : 
                j-=1
            
        elif fold[k]=="l":
            if flag == "right":
                j+=1
                flag = "up"
            elif flag=="left":
                j-=1
                flag = "down"
            elif flag == "up":
                i-=1
                flag="left"
            else : 
                i+=1
                flag = "right"
            
        else : #fold[k]=="r"
            if flag == "right":
                j-=1
                flag = "down"
            elif flag=="left":
                j+=1
                flag = "up"
            elif flag == "up":
                i+=1
                flag = "right"
            else : 
                i-=1
                flag = "left"
            
        # If it is an "h", we add it to the dictionnary because it can be used for the score
        if seq[k+1]=="h" or seq[k+1]=="H":
            if seq[k]=='h'or seq[k]=='H':
                grid[(i,j)]= k+1
            else : 
                grid[(i,j)]= k+1
        
     # Calculate the score : 
    keys = list(grid.keys())
    score = 0
    for k in range(len(keys)):
        i,j = keys[k]
        if (i-1,j) in grid and grid[i-1,j]%2!=grid[i,j]%2 and grid[i-1,j]!=grid[i,j]-1 and grid[i-1,j]!=grid[i,j]+1: 
            # the 2 last conditions are to avoid scoring H-h that are linked
            score+=1
        if (i,j+1) in grid and grid[i,j]%2 != grid[i,j+1]%2 and grid[i,j]!=grid[i,j+1]-1 and grid[i,j]!=grid[i,j+1]+1: 
            # the 2 last conditions are to avoid scoring H-h that are linked
            score +=1
    return(score)


#####################################################################
# Main
#####################################################################
seq_in = sys.argv[1]

#assign the odd/even labels
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
            left_part = (between // 2) - 1
            j = left_part
            directions = directions + "l"
            while(j > 0):
                directions = directions + "f"
                j = j - 1
            directions = directions + "rr"
            j = left_part
            while(j > 0):
                directions = directions + "f"
                j = j - 1
            directions = directions + "l"

#creating the main fold:
between_s1_s2 = pairs[len(pairs)-1][1] - pairs[len(pairs)-1][0] - 1
j = between_s1_s2 // 2
while(j > 0):
    directions = directions + "f"
    j = j-1
directions = directions + "rr"
j = (between_s1_s2 // 2) - 1 
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
            left_part = (between // 2) - 1
            j = left_part
            directions = directions + "l"
            while(j > 0):
                directions = directions + "f"
                j = j - 1
            directions = directions + "rr"
            j = left_part
            while(j > 0):
                directions = directions + "f"
                j = j - 1
            directions = directions + "l"

last_matched = pairs[0][1]


while last_matched != len(seq_in) - 1:
    directions = directions + "f"
    last_matched = last_matched + 1

#not a nice way how to fix one bug, would be nice to improve
if len(directions) == len(seq_in):
    directions = directions[0:len(directions) - 1]

#print(seq_in)
print("Fold:", directions)

score = score(seq_in, directions)
print("Score: ", score)

             











