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


##################Christian's code#######################
class HPFold:

    def __init__ (self, s):
        legal = {'h':'h', 'p':'p', 'H':'h', 'P':'p'}
        self.seq = []
        i = 1
        for c in s:
            if c in legal.keys():
                if legal[c] == 'h' and i % 2 == 0:
                    self.seq.append('H')
                else:
                    self.seq.append(legal[c])
                i = i + 1

    def __len__ (self):
        return len(self.seq)
        
    def SetRelFold (self, relfold):
        """
        Fold seq according to a description in relavtive format, i.e.
        a sequence of {f,l,r}'s which describe each step as (f)orward,
        (l)eft, or (r)ight.
        """
        turn = {'f':0, 'l':-1, 'r':1}
        direction = {0:'e', 1:'s', 2:'w', 3:'n'}
        absfold = []
        curr = 0
        for relstep in relfold:
            absstep = (curr + turn[relstep]) % 4
            absfold.append(direction[absstep])
            curr = absstep
        return self.SetAbsFold(absfold)

    def SetAbsFold (self, absfold):
        """
        Fold seq according to a description in absolute format, i.e.
        s sequence of {n,s,e,w}'s which describe each step as (n)orth,
        (s)outh, (e)ast, or (w)est.
        """
        self.legal_fold = (True, 0)
        self.grid = {}
        self.grid[0,0] = [0]
        i = j = self.min_i= self.max_i = self.min_j = self.max_j = 0
        k = 1
        for step in absfold:
            if step == 'n':
                i = i - 1
            elif step == 's':
                i = i + 1
            elif step == 'e':
                j = j + 1
            elif step == 'w':
                j = j - 1
            if (i,j) in self.grid.keys():
                self.legal_fold = (False, k)
                self.grid[i,j].append(k)
            else:
                self.grid[i,j] = [k]
            k = k + 1
            self.min_i = min(i, self.min_i)
            self.max_i = max(i, self.max_i)
            self.min_j = min(j, self.min_j)
            self.max_j = max(j, self.max_j)                        
        return self.legal_fold[0]

    def ContainNeighbors(self, l1, l2):
        """
        Returns true if there exists k1 in l1 and k2 in l2 such that
        abs(k1-k2) is 1, i.e. if the indices in l1 and l2 contain a
        pair of neighbors in seq.
        """
        res = False
        for k1 in l1:
            for k2 in l2:
                if abs(k1-k2) == 1:
                    res = True
        return res

    def ContainHHs(self, l1, l2):
        """
        Returns true if there exists k1 in l1 and k2 in l2 where there
        is a 'h' at position k1 and k2 in seq, i.e. if the indices in
        l1 and l2 contain a pair which can make a h-h bond.
        """
        res = False
        for k1 in l1:
            for k2 in l2:
                if (self.seq[k1] == "h" or self.seq[k1] == "H") and (self.seq[k2] == "h" or self.seq[k2] == "H"):
                    res = True
        return res

    def PrintFold (self):
        """
        Print fold and output its score
        """
        score = 0
        #print()
        for i in range(self.min_i, self.max_i+1):
            for j in range(self.min_j, self.max_j+1):
                if (i,j) in self.grid.keys():
                    l1 = self.grid[i,j]
                    if len(l1) == 1:
                        pass
                        #print(self.seq[l1[0]], end="")
                    else:
                        pass
                        #print("X", end="")
                    if (i,j+1) in self.grid.keys():
                        l2 = self.grid[i,j+1]
                        if self.ContainNeighbors(l1,l2):
                            pass
                            #print("-", end="")
                        elif self.ContainHHs(l1, l2):
                            #print("*", end="")
                            score = score + 1
                        else:
                            pass
                            #print(" ", end="")
                    else:
                        pass
                        #print(" ", end="")
                else:
                    pass
                    #print(".", end="")
                    #print(" ", end="")
            

            for j in range(self.min_j, self.max_j+1):
                if (i,j) in self.grid.keys() and (i+1,j) in self.grid.keys():
                    l1 = self.grid[i,j]
                    l2 = self.grid[i+1,j]
                    if self.ContainNeighbors(l1,l2):
                        pass
                        #print("|", end="")
                    elif self.ContainHHs(l1,l2):
                        #print("*", end="")
                        score = score + 1
                    else:
                        pass
                        #print(" ", end="")
                else:
                    pass
                    #print(" ", end="")
                #print(" ", end="")
            #print()

        if self.legal_fold[0]:
            print("Score: %d" % (score))
        else:
            print("Illegal fold after %d steps" % (self.legal_fold[1]))
       

#####################################################################
# Main
#####################################################################
seq_in = sys.argv[1]

#seq_in = "hphhhhhh"
#seq_in = "hhppppphhppphppphp"
#seq_in = "hphphhhppphhhhpphh"

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

seq = HPFold(seq_in)

seq.SetRelFold(directions)

seq.PrintFold()

             











