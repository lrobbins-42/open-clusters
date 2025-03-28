# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 13:09:09 2025

@author: lucie
"""

import numpy as np
from iminuit import Minuit
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
#import iminuit
#from iminuit import Minuit
import read_mist_models

# Define chi-squared function
def chi2(parin):
    mu = parin[0]
    Av = parin[1]
    
    colour_excess = Av/Rv
    # Transform isochrone
    B_V_transformed = BminV_iso + colour_excess
    V_transformed = V_iso + mu + Av

    # Compute squared residuals by finding closest isochrone point to each observation
    chi_sq = 0
    for B_Vo, Vo in zip(B_min_V_cal, V_calibrated):
        distances = np.sqrt((B_V_transformed - B_Vo)**2 + (V_transformed - Vo)**2)
        
        # Find the minimum distance (closest point on the isochrone) and add squared residual
        idx_min = np.argmin(distances)
        chi_sq += distances[idx_min]**2

    return chi_sq

folder_name = 'FINAL/NGC 744'
cluster_name = 'NGC 744'
name = 'NGC744'
file_path = '/Users/lucie/OneDrive/Documents/Uni Y4/Major Project/Calibration/NEW/'+folder_name+'/'  
 
excel_file = file_path+name+'_calibrated_magnitudes.xlsx'
#excel_file = file_path+name+'_no_outliers.xlsx'
M1 = pd.read_excel(excel_file, sheet_name=0)  #READ B
M2 = pd.read_excel(excel_file, sheet_name=1)  #READ V

#read in NGC 744 parameters to calculate calibrated magnitudes
param_file = '/Users/lucie/OneDrive/Documents/Uni Y4/Major Project/Calibration/NEW/FINAL/NGC 744/NGC744_calibration_parameters.xlsx'
P = pd.read_excel(param_file, sheet_name=0)
MLE = list(P['MLE'])
MLE_Error = list(P['MLE Error'])

Alpha = MLE[0]
Beta = MLE[1]
b_zero = MLE[2]
v_zero = MLE[3]

alpha_error = MLE_Error[0]
beta_error = MLE_Error[1]
b0_error = MLE_Error[2]
v0_error = MLE_Error[3]

zero_point = np.array([[b_zero], [v_zero]])

B_calibrated = []
V_calibrated = []
B_min_V_inst = []
B_min_V_cal = []
B_min_V_SIM = []
B_min_V_errors = []
V_cal_errors = []

B_inst = list(M1['B inst'])
B_error = list(M1['B error'])
B_sim = list(M1['B cal'])
V_inst = list(M2['V inst'])
V_error = list(M2['V error'])
V_sim = list(M2['V cal'])

for b_i,v_i in enumerate(V_inst): #find calibrated BVs based on the inst. mags and 
                                  #the parameters we've estimated
    c = np.array([[1-Alpha, Alpha], [Beta, 1-Beta]])
    d = np.vstack([B_inst[b_i], v_i])
    x = np.matmul(c,d)+zero_point
    #solves for b, v as in notes eq.1
    b_cal = x[0] 
    v_cal = x[1]
    B_calibrated.append(b_cal)
    V_calibrated.append(v_cal)
    
    #find SIMBAD b-v values
    sim_b_min_v = B_sim[b_i]-V_sim[b_i]
    B_min_V_SIM.append(sim_b_min_v)
    
    #create calibrated b-v array for CMD
    min_item = B_inst[b_i]-v_i
    B_min_V_inst.append(min_item)
    min_cal_item = ((1-Alpha-Beta)*(min_item))+(b_zero-v_zero)
    B_min_V_cal.append(min_cal_item)
    
    #find errors on each v_cal point
    sigma_v_cal = ((Beta**2)*B_error[b_i]**2)+(((1-Beta)**2)*V_error[b_i]**2)
    V_cal_errors.append(np.sqrt(sigma_v_cal))
    
    #find errors on each B-V point
    #b_i is just the index so can be used for both here
    B_min_V_error = np.sqrt(B_error[b_i]**2 + V_error[b_i]**2)    
    sigma_b_min_v = (1-Alpha-Beta)*B_min_V_error
    B_min_V_errors.append(sigma_b_min_v)

V_idxs = list(M2['Star Index'])
#plot the CMD
figA, (ax1) = plt.subplots(1,1, figsize=(5,6))
ax1.errorbar(B_min_V_cal, V_calibrated, yerr=V_cal_errors, xerr=B_min_V_errors, fmt='o', color='b', ecolor='black', mec='black', label='Observational Data')
ax1.scatter(B_min_V_SIM, V_sim, color='r',label='SIMBAD Data')
ax1.set_ylabel(r'$v$')
ax1.set_xlabel(r'$(b-v)$')
ax1.set_title(cluster_name+' CMD')
ax1.invert_yaxis()
for i, txt in enumerate(V_idxs):
    plt.annotate(txt, (B_min_V_cal[i], V_calibrated[i]), textcoords="offset points",
                 xytext=(5,5), ha='right',  fontsize=6)
ax1.legend()#loc='upper left'
    
#plot B-V in each axis to test if straight line
#plot inst against the SIMBAD magnitudes
figB, (ax1) = plt.subplots(1,1, figsize=(5,6))
ax1.scatter(B_min_V_SIM, B_min_V_inst, s=12)
#create straight line fit
xmin = min(B_min_V_SIM)
xmax = max(B_min_V_SIM)
x_points = np.linspace(xmin,xmax,100)
y_plot = (x_points-(b_zero-v_zero))/(1-Alpha-Beta)
#ax1.plot(x_points, y_plot+0.05, color="r")
ax1.plot(x_points, y_plot, color="r", label=r'Ideal $(B-V)$')
ax1.set_ylabel(r'$(B-V)$')
ax1.set_xlabel(r'$(b-v)_{\mathrm{SIMBAD}}$')
ax1.set_title(cluster_name+' B-V Colour Index')
ax1.legend(loc='upper left', bbox_to_anchor=(0.05, 0.95))

#ISOCHRONE FITTING BELOW

#get isochrone magnitude data from MIST
file3 = read_mist_models.ISOCMD('MIST_iso_89.iso.cmd') #NGC 744, M34,5,7
file4 = read_mist_models.ISOCMD('MIST_iso_78.iso.cmd') #M36

n_decimal = 3
isocmd = file3
n = [round(x,2) for x in isocmd.ages]
n_short = [n[n_decimal]] #M34,5,7, NGC 744
#n_short = [n[n_decimal]] #M36 for file4

print('version: ', isocmd.version)
print('photometric system: ', isocmd.photo_sys)
print('abundances: ', isocmd.abun)
Z_solar = isocmd.abun['Zinit']
print('rotation: ', isocmd.rot)
print('ages: ', n_short)
print('number of ages: ', isocmd.num_ages)
print('available columns: ', isocmd.hdr_list)
print('Av extinction: ', isocmd.Av_extinction)

for item in n_short:
    age_ind = isocmd.age_index(item) #returns the index for the desired age
    B = isocmd.isocmds[age_ind]['Bessell_B'] #Bessell filters
    V = isocmd.isocmds[age_ind]['Bessell_V']
    V = np.array(V)
    BminV = np.array(B-V)

V_iso = []
BminV_iso = []

for i,v in enumerate(V):    
    if (v >= -3 and v <= 8) and (BminV[i] >= -0.2 and BminV[i] <= 1.8):
    #if (v >= -9 and v <= 8) and (BminV[i] >= -2 and BminV[i] <= 2):
        V_iso.append(v)
        BminV_iso.append(BminV[i])

# Initial parameter guesses
mu_init = 10  # Approximate distance modulus 10 for NGC 744, 8 for rest
Av_init = 0.8       # Estimated extinction in V
Rv = 3.1

parin = np.array([mu_init, Av_init])
parname = ['mu','Av']
parstep = np.array([0.1, 0.1]) #GUESSED BY ORDER OF MAG.
parfix  = [False, False] 
parlim = [(7, 12), (0, 3)] 
#parlim = [(7, 11), (None, 2)] 
# Use iminuit to minimize chi-squared
m = Minuit(chi2, parin, name=parname)

m.errors = parstep
m.limits = parlim
m.fixed = parfix
# Minimize
m.migrad()

# Get best-fit values and uncertainties
mu_fit, Av_fit = m.values["mu"], m.values["Av"]
mu_err, Av_err = m.errors["mu"], m.errors["Av"]

print(f"Best-fit parameters:")
print(f"  mu = {mu_fit:.3f} ± {mu_err:.3f}")
print(f"  Av = {Av_fit:.3f} ± {Av_err:.3f}")
print(f"  E(B-V) = {Av_fit/Rv:.3f} ± {Av_err/Rv:.3f}")

V = np.array(V_iso) +mu_fit +Av_fit
B_V = np.array(BminV_iso) +(Av_fit/Rv)

fig, (ax1) = plt.subplots(1,1, figsize=(5,6))
ax1.plot(B_V, V, label=fr'$n$={n_short[0]}', color='r') #plot isochrone
ax1.errorbar(B_min_V_cal, V_calibrated, yerr=V_cal_errors, xerr=B_min_V_errors, fmt='o', color='b', ecolor='black', mec='black', label='Observational \nData')#with star \n index from \n V image
ax1.set_ylabel(r'$v$')
ax1.set_xlabel(r'$(b-v)$')
ax1.set_title(cluster_name+' CMD with MIST Isochrones')
ax1.invert_yaxis()
for i, txt in enumerate(V_idxs):
    plt.annotate(txt, (B_min_V_cal[i], V_calibrated[i]), textcoords="offset points",
                 xytext=(5,5), ha='right',  fontsize=6)
fig.legend(loc='lower right', bbox_to_anchor=(0.9, 0.1))
fig.show()
#plt.pause(0.1)

print('\n'+'Cluster Parameters')
distance_estimate = 10**(1+(mu_fit/5))
print(f'Distance to cluster: {distance_estimate:.2f} Pc')
print(f'Z: {Z_solar}')
age_estimate = (10**n_short[0])/1E6 #in million yrs
print(f'Age of cluster: {age_estimate:.2f} Million Years \n')

#write info to text file
#output_file = file_path+name+f'_cluster_values_v{n_decimal}.txt'
#output_file = file_path+name+f'_cluster_values_v{n_decimal}_NO.txt'
#with open(output_file, 'w') as f:
 #   f.write(f'Distance to cluster: {distance_estimate:.2f} Pc \n')
  #  f.write(f'Distance modulus (mu): {mu_fit:.3f} ± {mu_err:.3f} \n')
   # f.write(f'Av: {Av_fit:.3f} ± {Av_err:.3f} \n')
    #f.write(f'E(B-V): E(B-V) = {Av_fit/Rv:.3f} ± {Av_err/Rv:.3f} \n')
    #f.write(f'Z: {Z_solar} \n')
    #f.write(f'Age of cluster: {age_estimate:.2f} Million Years')
#print(f"Cluster values have been written to {output_file}")