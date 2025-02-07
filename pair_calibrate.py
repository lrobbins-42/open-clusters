# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 14:39:45 2024

@author: lucie
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

cluster_name = 'NGC 744'
name = 'NGC744' #change based on image
file_path = '/Users/lucie/OneDrive/Documents/Uni Y4/Major Project/Calibration/NEW/C-'+cluster_name+'/'

#both simbad files are IDENTICAL so just use one for both B,V
S = pd.read_excel(file_path+name+"B_simbad_magnitudes.xlsx") #SIMBAD
D = pd.read_excel(file_path+name+"_dao_magnitudes_B.xlsx") #DAO B mag

S_RA = list(S['RA (deg)'])
S_DEC = list(S['DEC (deg)']) 
D_RA = list(D['RA'])
D_DEC = list(D['DEC'])
threshold = 0.005 

B_simbad_indices = []
B_dao_indices = [] #list of which indexed star is star 1,2,3.. in SIMBAD

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
                B_dao_indices.append(i)
                B_simbad_indices.append(closest_match)

print('Threshold = '+str(threshold))
print("B DAO Indices:", B_dao_indices)
print("B SIMBAD Indices:", B_simbad_indices)

#REPEAT FOR V MAGNITUDES
D_V = pd.read_excel(file_path+name+"_dao_magnitudes_V.xlsx") #DAO V mag
D_V_RA = list(D_V['RA'])
D_V_DEC = list(D_V['DEC']) 
S_V_RA = list(S['RA (deg)'])
S_V_DEC = list(S['DEC (deg)'])    

V_simbad_indices = []
V_dao_indices = [] #list of which indexed star is star 1,2,3.. in SIMBAD

for i, RA_value in enumerate(D_V_RA): #testing the DAO RA against the simbad RA (r)
    RAmatches = [(j, r) for j, r in enumerate(S_V_RA) if abs(RA_value - r) <= threshold]
    if RAmatches:
        RAmatch_indices = [match[0] for match in RAmatches]
        DEC_value = D_V_DEC[i] #testing the DAO DEC against the simbad DEC (d)
        DECmatches = [(j, d) for j, d in enumerate(S_V_DEC) if abs(DEC_value - d) <= threshold]
        if DECmatches:
            DECmatch_indices = [match[0] for match in DECmatches]
            common_indices = set(RAmatch_indices).intersection(DECmatch_indices)
            if common_indices:
                # Find the closest match based on RA and DEC combined
                closest_match = min(
                    common_indices,
                    key=lambda j: abs(RA_value - S_V_RA[j]) + abs(DEC_value - S_V_DEC[j]),
                )
                # Append the closest match only
                V_dao_indices.append(i)
                V_simbad_indices.append(closest_match)
  
print("V DAO Indices:", V_dao_indices)
print("V SIMBAD Indices:", V_simbad_indices)

#NOW COMPARE THE ARRAYS TO FIND THE MATCHED PAIRS OF DAO STARS AND WHICH SIMBAD
#STARS THEY CORRESPOND TO. HAS TO BE SEPARATE FOR B AND V AS THE SIMBAD INDICES
#ARE DIFFERENT FOR THE B AND V IF YOU LOOK AT IT CAREFULLY
simbad_matched_pairs = []
sim_idxs = []
dao_B_pairs = []
dao_V_pairs = []
for s_idx, s1 in enumerate(V_simbad_indices):
    if s1 in V_simbad_indices and s1 in B_simbad_indices:
        simbad_matched_pairs.append(s1)  
        sim_idxs.append(s_idx) #corresponding indices for SIMBAD array search
for s2 in sim_idxs:
    dao_B_pairs.append(B_dao_indices[s2])
    dao_V_pairs.append(V_dao_indices[s2])

print('\n'+'Below are the star indices that we want for the pair calibration')     
print('SIMBAD Matches', simbad_matched_pairs)
print('DAOPHOT B Matches', dao_B_pairs)
print('DAOPHOT V Matches', dao_V_pairs)

DAO_B_mag = list(D['B']) #reading the magnitudes from the excel sheets
SIM_B_mag = list(S['B'])
DAO_V_mag = list(D_V['V'])
SIM_V_mag = list(S['V'])
B_err = list(D['B Error']) #errors on magnitudes from DAO
V_err = list(D_V['V Error'])

B_inst = [] #instrumental DAOPHOT magnitudes
V_inst = [] 
B_errors = [] #and corresponding errors
V_errors = []
B_cal = []  #calibrated SIMBAD magnitudes for B
V_cal = [] #and V

for n in dao_B_pairs:
    print(n)
    B_inst.append(DAO_B_mag[n])  
    B_errors.append(B_err[n])
    
for p in dao_V_pairs:
    V_inst.append(DAO_V_mag[p]) 
    V_errors.append(V_err[p])
    
for m in simbad_matched_pairs:
    B_cal.append(SIM_B_mag[m])
    V_cal.append(SIM_V_mag[m])
      
print('V inst:', V_inst)
print('V cal:', V_cal)
print('B inst:', B_inst)
print('B cal:', B_cal)
    

#WRITE INFO TO FILE
B_data = {
    'B inst': B_inst,
    'B error': B_errors,
    'B cal': B_cal
}
#treat B and V separately as arrays may be different lengths
V_data = {
    'V inst': V_inst,
    'V error': V_errors,
    'V cal': V_cal
}


dfB = pd.DataFrame(B_data)
dfV = pd.DataFrame(V_data)
excel_file = file_path+name+'_calibrated_magnitudes.xlsx'
with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
    dfB.to_excel(writer, index=False, sheet_name='Sheet 1', float_format='%.4f')
    dfV.to_excel(writer, index=False, sheet_name='Sheet 2', float_format='%.4f')

    worksheet1 = writer.sheets['Sheet 1']
    worksheet2 = writer.sheets['Sheet 2']
  
    for col_num, value in enumerate(dfB.columns.values):
        worksheet1.column_dimensions[chr(65 + col_num)].width = 12
        
    for col_num, value in enumerate(dfB.columns.values):
        worksheet2.column_dimensions[chr(65 + col_num)].width = 12
    print(f'\nSummary data has been written to {excel_file}')

#PLOTTING INST-CAL 
plt.ion()   
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(5,6))
ax1.scatter(V_inst, V_cal, s=12)
A,B = np.polyfit(V_inst, V_cal, 1)
best_fit_line = (A*np.array(V_inst)) + B
ax1.plot(V_inst, best_fit_line, color="r")
ax1.set_xlabel(r'$V_{inst}$')
ax1.set_ylabel(r'$V_{cal}$')
ax1.set_title(cluster_name+' Star Pairs for V,B Magnitudes')

plt.subplots_adjust(hspace=0.3)

ax2.scatter(B_inst, B_cal, s=12)
C,D = np.polyfit(B_inst, B_cal, 1)
best_fit_line = (C*np.array(B_inst)) + D
ax2.plot(B_inst, best_fit_line, color="r")
ax2.set_xlabel(r'$B_{inst}$')
ax2.set_ylabel(r'$B_{cal}$')

#PLOTTING B-V INSTRUMENTAL and SIMBAD CMD
BV_i = np.array(B_inst)-np.array(V_inst)
BV_s = np.array(B_cal)-np.array(V_cal)
fig2, (ax1, ax2) = plt.subplots(2,1, figsize=(5,6))
ax1.plot(BV_i, V_inst, '-b+')
ax1.set_xlabel(r'$V_{inst}$')
ax1.set_ylabel(r'$B_{inst}-V_{inst}$')
ax1.set_title(cluster_name+' Trial CMDs')
ax1.invert_yaxis()

plt.subplots_adjust(hspace=0.3)

ax2.plot(BV_s, V_cal, '-b+')
ax2.set_xlabel(r'$V_{cal}$')
ax2.set_ylabel(r'$B_{cal}-V_{cal}$')
ax2.invert_yaxis()