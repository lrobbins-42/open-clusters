# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:39:45 2024

@author: lucie
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

S = pd.read_excel("simbad_magnitudes.xlsx") #SIMBAD
D = pd.read_excel("dao_magnitudes_v.xlsx") #DAO

S_RA = list(S['RA (deg)'])
S_DEC = list(S['DEC (deg)'])
D_RA = list(D['RA'])
D_DEC = list(D['DEC'])
threshold = 0.2 #0.02 is BEST VALUE SO FAR

selected_values = []
simbad_indices = []
dao_indices = [] #list of which indexed star is star 1,2,3.. in SIMBAD

for i, RA_value in enumerate(D_RA): #testing the DAO RA against the simbad RA (r)
    RAmatches = [(j, r) for j, r in enumerate(S_RA) if abs(RA_value - r) <= threshold]
    if RAmatches:
        RAmatch_indices = [match[0] for match in RAmatches]
        DEC_value = D_DEC[i] #testing the DAO DEC against the simbad DEC (d)
        DECmatches = [(j, d) for j, d in enumerate(S_DEC) if abs(DEC_value - d) <= threshold]
        if DECmatches:
            DECmatch_indices = [match[0] for match in DECmatches]
            common_indices = set(RAmatch_indices).intersection(DECmatch_indices)
            if common_indices:
                # Find the closest match based on RA and DEC combined
                closest_match = min(
                    common_indices,
                    key=lambda j: abs(RA_value - S_RA[j]) + abs(DEC_value - S_DEC[j]),
                )
                # Append the closest match only
                dao_indices.append(i)
                simbad_indices.append(closest_match)

print('Threshold = '+str(threshold))
print("DAO Indices:", dao_indices)
print("SIMBAD Indices:", simbad_indices)

            
