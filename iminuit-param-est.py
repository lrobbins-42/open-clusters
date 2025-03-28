# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:51:23 2025

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
    x = np.matmul(a_inv,y)
    B = x[0]
    V = x[1]

    return B, V

def chi_squared(parin):  
    alpha = parin[0]
    beta = parin[1]
    b0 = parin[2]
    v0 = parin[3]

    MBi = np.array(M1['B inst'])
    MBc = np.array(M1['B cal']) 
    #sigma_B = 50*np.array(M1['B error'])
    sigma_B = 1*np.array(M1['B error']) 
    
    MVi = np.array(M2['V inst'])
    MVc = np.array(M2['V cal'])
    #sigma_V = 50*np.array(M2['V error'])
    sigma_V = 1*np.array(M2['V error'])

    b_terms = []
    v_terms = []    
    N = len(MBi)
    #print('ndof: ',(2*N)-4)
    #sigma_B = np.full(N, 0.6)
    #sigma_V = np.full(N, 0.6)
    
    for x in range(N):
        B, V = BV_matrix(MBc[x], MVc[x], alpha, beta, b0, v0) #solve matrix
        #print('B inst: ',B)
        #print('V inst: ',V)
        b_terms.append(((MBi[x]-B)**2)/(sigma_B[x]**2))
        v_terms.append(((MVi[x]-V)**2)/(sigma_V[x]**2))
    
    B_term = np.sum(b_terms) #sum all ith terms
    V_term = np.sum(v_terms)
    #print('B: ',B_term)
    #print('V: ',V_term)
    chi2 = B_term + V_term #final chi^2
    #print('chi2: ',chi2)
    
    return chi2

folder_name = 'FINAL/NGC 744'
cluster_name = 'NGC 744'
name = 'NGC744'
file_path = '/Users/lucie/OneDrive/Documents/Uni Y4/Major Project/Calibration/NEW/'+folder_name+'/'  
 
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
#parfix  = [False, False, True, True]  
#parfix  = [True, True, False, False]  
parlim = [(-10, 10), (-10, 10), (None, None), (None, None)] 
#parlim = [(-0.5, 0.5), (-0.5, 0.5), (None, None), (None, None)] 
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
        f.write(f"{parname[n]}: {MLE[n]} Error: {sigmaMLE[n]} \n")

print(f"Parameter values have been written to {output_file}")

param_data = {
    'MLE': MLE,
    'MLE Error': sigmaMLE
}
df = pd.DataFrame(param_data)
excel_file = file_path+name+'_calibration_parameters.xlsx'
with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
    df.to_excel(writer, index=False, sheet_name='Sheet 1', float_format='%.8f')
    worksheet1 = writer.sheets['Sheet 1']
    for col_num, value in enumerate(df.columns.values):
        worksheet1.column_dimensions[chr(65 + col_num)].width = 12
    print(f'\nSummary data has been written to {excel_file}')

