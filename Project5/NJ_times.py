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
    order = "python3 project5_numpy.py " + file
    start_time = time.time()
    os.system(order)
    times.append(time.time() - start_time)
    print(i,times[i])
print (times)

# times = [0.6804704666137695, 311.48744106292725, 123.53081893920898, 58.615458965301514, 32.855236530303955, 16.757259368896484, 7.541655778884888, 2.7409214973449707, 1686.475982427597, 1435.3768684864044, 1271.1676616668701, 1000.0718369483948, 872.0920331478119, 639.3578417301178]