### time comparisons RapidNJ project 5

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
    order = "./rapidnj.exe " + file
    start_time = time.time()
    os.system(order)
    times.append(time.time() - start_time)
print (times)

# [0.765866756439209, 0.8584439754486084, 0.8805592060089111, 1.065643548965454, 1.2064378261566162, 1.2218735218048096, 0.056989431381225586, 0.07906603813171387, 0.11399245262145996, 0.13755559921264648, 0.1995706558227539, 0.2588193416595459, 0.35033106803894043, 0.030724525451660156]