### time comparisons project 5

import os
import time

dist_matrix = ['1347_FAINT.phy',
 '1493_Fe-ADH.phy',
 '1560_Ferritin.phy',
 '1689_FGGY_N.phy',
 '1756_FAD_binding_3.phy',
 '1849_FG-GAP.phy',
 '214_Arena_glycoprot.phy',
 '304_A1_Propeptide.phy',
 '401_DDE.phy',
 '494_Astro_capsid.phy',
 '608_Gemini_AL2.phy',
 '777_Gemini_V1.phy',
 '877_Glu_synthase.phy',
 '89_Adeno_E3_CR1.phy']

times = []
for i in range(len(dist_matrix)):
    file = dist_matrix[i]
    order = "python3 NJ.py " + file
    start_time = time.time()
    os.system(order)
    times.append(time.time() - start_time)
print (times)