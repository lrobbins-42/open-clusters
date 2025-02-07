# -*- coding: utf-8 -*-
"""
Created on Sat Feb  1 16:04:15 2025

@author: lucie
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import iminuit
from iminuit import Minuit
import read_mist_models

print("iminuit version:", iminuit.__version__)  # need 2.x

#function to solve the BV matrix in notes eq.2
def BV_matrix(MBc, MVc, alpha, beta, b0, v0):  
    a = np.array([[1-alpha, alpha], [beta, 1-beta]])
    y = np.vstack([MBc-b0, MVc-v0])
    a_inv = np.linalg.inv(a) #1/det * inverted matrix elements
    x = np.matmul(a,y)
    B = x[0]
    V = x[1]

    return B, V

def chi_squared(parin):  
    alpha = parin[0]
    beta = parin[1]
    b0 = parin[2]
    v0 = parin[3]
    
    MBi = list(M1['B inst'])
    MBc = list(M1['B cal']) 
    sigma_B = list(M1['B error'])
    
    MVi = list(M2['V inst'])
    MVc = list(M2['V cal'])
    sigma_V = list(M2['V error'])

    b_terms = []
    v_terms = []    
    N = len(MBi)
    
    for x in range(N):
        B, V = BV_matrix(MBc[x], MVc[x], alpha, beta, b0, v0) #solve matrix
        b_terms.append(((MBi[x]-B)**2)/(sigma_B[x]**2))
        v_terms.append(((MVi[x]-V)**2)/(sigma_V[x]**2))
    
    B_term = np.sum(b_terms) #sum all ith terms
    V_term = np.sum(v_terms)
    chi2 = B_term + V_term #final chi^2
    
    return chi2

cluster_name = 'NGC 744' #CHANGE THESE FOR EVERY IMAGE
name = 'NGC744'
file_path = '/Users/lucie/OneDrive/Documents/Uni Y4/Major Project/Calibration/NEW/C-'+cluster_name+'/'  
 
excel_file = file_path+name+'_calibrated_magnitudes.xlsx'
M1 = pd.read_excel(excel_file, sheet_name=0)  #READ B
M2 = pd.read_excel(excel_file, sheet_name=1)  #READ V
  
#set initial guesses 
alpha = 0
beta = 0
b0 = 18
v0 = 18

#iminuit stuff
parin = np.array([alpha, beta, b0, v0])
parname = ['alpha','beta','b0','v0']
parstep = np.array([0.1, 0.1, 1., 1.]) #GUESSED BY ORDER OF MAG.
parfix  = [False, False, False, False]  
parlim = [(-10, 10), (-10, 10), (None, None), (None, None)] 
m = Minuit(chi_squared, parin, name=parname)

m.errors = parstep
m.limits = parlim
m.fixed = parfix
m.migrad()                                        #minimised chi2
MLE = m.values                                    #max-likelihood estimates
sigmaMLE = m.errors                               #standard deviations
cov = m.covariance                                #covariance matrix
rho = m.covariance.correlation()  

#output results to text file
print(f"Maximum Likelihood Estimates: {MLE} \n")
print(f"Parameter Errors: {sigmaMLE} \n")
print(f"Covariance Matrix: {cov}")
print(f"Correlation Matrix: {rho}")

output_file = file_path+name+'_parameter_values.txt'
with open(output_file, 'w') as f:
    #write the parameter names and their corresponding values
    for n in range(4):
        f.write(f"{parname[n]}: {MLE[n]} \n")

print(f"Parameter values have been written to {output_file}")

#use parameters to calculate calibrated magnitudes
Alpha = MLE[0]
Beta = MLE[1]
b_zero = MLE[2]
v_zero = MLE[3]
zero_point = np.vstack([b_zero, v_zero])

B_calibrated = []
V_calibrated = []
B_min_V_inst = []
B_min_V_cal = []

B_inst = list(M1['B inst'])
V_inst = list(M2['V inst'])
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
    
    #create calibrated B-V array for CMD
    min_item = B_inst[b_i]-v_i
    B_min_V_inst.append(min_item)
    min_cal_item = ((1-(Alpha-Beta))*(min_item))+(b_zero-v_zero)
    B_min_V_cal.append(min_cal_item)

#plot calibrated BV
plt.ion()   
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(5,6))
ax1.scatter(V_inst, V_calibrated, s=12)
A1,B1 = np.polyfit(V_inst, V_calibrated, 1)
best_fit_line = (A1*np.array(V_inst)) + B1
ax1.plot(V_inst, best_fit_line, color="r")
ax1.set_xlabel(r'$V_{inst}$')
ax1.set_ylabel(r'$V_{cal}$')
ax1.set_title(cluster_name+' Calibrated V,B Magnitudes')

plt.subplots_adjust(hspace=0.3)

ax2.scatter(B_inst, B_calibrated, s=12)
C1,D1 = np.polyfit(B_inst, B_calibrated, 1)
best_fit_line = (C1*np.array(B_inst)) + D1
ax2.plot(B_inst, best_fit_line, color="r")
ax2.set_xlabel(r'$B_{inst}$')
ax2.set_ylabel(r'$B_{cal}$')

#plot B-V in each axis to test if straight line
figB, (ax1) = plt.subplots(1,1, figsize=(5,6))
ax1.scatter(B_min_V_inst, B_min_V_cal, s=12)
ax1.set_xlabel(r'$(B-V)_{inst}$')
ax1.set_ylabel(r'$(B-V)_{cal}$')
ax1.set_title(cluster_name+' Calibrated B-V Colour Index')


