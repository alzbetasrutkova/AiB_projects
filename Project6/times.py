# Project 4

import os
import time

list_seqs= ["hhppppphhppphppphp",
            "hphphhhppphhhhpphh",
            "phpphphhhphhphhhhh",
            "hphpphhphpphphhpphph",
            "hhhpphphphpphphphpph",
            "hhpphpphpphpphpphpphpphh",
            "pphpphhpppphhpppphhpppphh",
            "ppphhpphhppppphhhhhhhpphhpppphhpphpp",
            "pphpphhpphhppppphhhhhhhhhhpppppphhpphhpphpphhhhh",
            "hhphphphphhhhphppphppphpppphppphppphphhhhphphphphh",
            "pphhhphhhhhhhhppphhhhhhhhhhphppphhhhhhhhhhhhpppphhhhhhphhphp",
            "hhhhhhhhhhhhphphpphhpphhpphpphhpphhpphpphhpphhpphphphhhhhhhhhhhh",
            "hhhhpppphhhhhhhhhhhhpppppphhhhhhhhhhhhppphhhhhhhhhhhhppphhhhhhhhhhhhppphpphhpphhpphph",
            "pppppphphhppppphhhphhhhhphhpppphhpphhphhhhhphhhhhhhhhhphhphhhhhhhppppppppppphhhhhhhpphphhhpppppphphh",
            "ppphhpphhhhpphhhphhphhphhhhpppppppphhhhhhpphhhhhhppppppppphphhphhhhhhhhhhhpphhhphhphpphphhhpppppphhh"
    ]

times = []
for i in range(len(list_seqs)):
    file = list_seqs[i]
    order = "python3 approx_fold.py " + file
    start_time = time.time()
    os.system(order)
    times.append(time.time() - start_time)
    print(i,times[i])
print (times)

