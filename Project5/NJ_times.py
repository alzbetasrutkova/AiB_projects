### time comparisons project 5

import os
import time

dist_matrix = ['89_Adeno_E3_CR1.phy',
 '877_Glu_synthase.phy',
 '777_Gemini_V1.phy',
 '608_Gemini_AL2.phy',
 '494_Astro_capsid.phy',
 '401_DDE.phy',
 '304_A1_Propeptide.phy',
 '214_Arena_glycoprot.phy',
 '1849_FG-GAP.phy',
 '1756_FAD_binding_3.phy',
 '1689_FGGY_N.phy',
 '1560_Ferritin.phy',
 '1493_Fe-ADH.phy',
 '1347_FAINT.phy']

times = []
for i in range(len(dist_matrix)):
    file = dist_matrix[i]
    order = "python3 NJ.py " + file
    start_time = time.time()
    os.system(order)
    times.append(time.time() - start_time)
    print(i,times[i])
print (times)