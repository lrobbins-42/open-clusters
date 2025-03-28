# -*- coding: utf-8 -*-
"""
Created on Sun Mar 16 23:26:06 2025

@author: lucie
"""
#handling both filters at once

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.time import Time
from scipy.optimize import curve_fit

def variable_analysis(N,period,FILTER,star_idxs,ref_idxs):
    M_diff = []
    M_error = []
    obs_times = []
    days = []
    
    for n in range(N):
        if star_idxs[n] >0 and ref_idxs[n] >0: 
            #print(f'{n} accepted')
            image_file = f'ngc744_{n}-0001_60s_{FILTER}.fits' 
            hdul = fits.open(file_path+image_file)
            header = hdul[0].header
            #get observation time
            date_time = header.get('DATE-OBS')
            hdul.close()
            #print(date_time)
            
            time_obs = Time(date_time, format="isot", scale="utc")
            day = time_obs.datetime.day
            days.append(day)
            #extract hours, minutes, seconds (ignoring milliseconds)
            hours = time_obs.datetime.hour
            minutes = time_obs.datetime.minute
            seconds = time_obs.datetime.second
            #convert time to total seconds
            total_seconds = (hours*3600) +(minutes*60) +seconds
            obs_times.append(total_seconds)
            
            #read inst magnitudes from the files
            excel_file = f'NGC744_{n}_dao_magnitudes_{FILTER}.xlsx' 
            M = pd.read_excel(file_path+excel_file, sheet_name=0)
            flux = np.array(M[f'{FILTER}'])
            error = np.array(M[f'{FILTER} Error'])
            
            star_mag = flux[star_idxs[n]]
            star_error = error[star_idxs[n]]
            ref_mag = flux[ref_idxs[n]]
            ref_error = error[ref_idxs[n]]
            error_prop = np.sqrt(star_error**2 + ref_error**2)
            
            M_diff.append(star_mag-ref_mag)
            M_error.append(error_prop)
    
    zero_time = [days[0],obs_times[0]] 
    times_since_zero = []
    
    for i,time in enumerate(obs_times):
        diff = days[i]-zero_time[0]
        day_secs = (diff*24*60*60)+time-zero_time[1] #seconds since the zero day
        times_since_zero.append(day_secs)
    
    M_diff = np.array(M_diff)
    TmodP = np.mod(times_since_zero,period)
    
    return M_diff, M_error, TmodP

# Improved model with sinusoidal eclipses + ellipsoidal variation + phase shift
def improved_eclipse_model(time_t,FILTER):  
    if FILTER == 'V':
        F_base = 0.17
        A1 = 0.09  #primary eclipse depth
        A2 = 0.09 #secondary eclipse depth
        A3 = 0.04 #third coeff guess
        a = c = 4
        b = a/2
        t_prim = 2500 
        t_sec = 1300 
        phi = 1800
    
    if FILTER == 'R':
        F_base = 0.022 
        A1 = 0.09  #primary eclipse depth
        A2 = 0.09 #secondary eclipse depth
        A3 = 0.04 #third coeff guess
        a = c = 4
        b = a/2
        t_prim = 2500 
        t_sec = 1300 
        phi = 2000
    
    primary_eclipse = A1* np.cos(a*np.pi * (time_t-t_prim + phi) / period)
    secondary_eclipse = A2 * np.cos(b*np.pi * (time_t-t_sec + phi) / period)
    third_term = A3 * np.cos(c*np.pi*(time_t + phi)/period)
    F_t = F_base -primary_eclipse + secondary_eclipse +third_term
    
    return F_t

folder_name = 'NGC744-Variability1'
cluster_name = 'NGC 744'
name = 'NGC744'
file_path = '/Users/lucie/OneDrive/Documents/Uni Y4/Major Project/Calibration/NEW/'+folder_name+'/'  

FILTER1 = 'V'
V1_idxs = np.array([64,0,65,63,62,63,0,49,53,48,0,41,39,53,44,0,0,46,50,50,35])-1 #index of star V1 in the V DAO images
V_ref_idxs = np.array([46,0,48,47,47,48,0,37,38,31,0,26,24,36,29,0,0,30,33,34,23])-1 #test ref star

R1_idxs = np.array([81,71,61,65,73,71,59,57,56,54,53,49,47,0,55,58,58,48,47,0,39])-1
R_ref_idxs = np.array([58,53,43,48,56,56,42,43,42,35,34,34,32,0,39,39,39,31,30,0,26])-1
FILTER2 = 'R'

N = 21
period = 0.3416*24*60*60 #0.3416 days periodicity of V1

#RUN V
V_diff, V_error, V_TmodP = variable_analysis(N,period,FILTER1,V1_idxs,V_ref_idxs)
#RUN R
R_diff, R_error, R_TmodP = variable_analysis(N,period,FILTER2,R1_idxs,R_ref_idxs)

# Generate fitted model curve
time_t = np.linspace(0-(0.3*period), 0+(0.7*period), 1000)
#time_t = np.linspace(0, 3*period, 1000)
F_t_V = improved_eclipse_model(time_t,FILTER1)
F_t_R = improved_eclipse_model(time_t,FILTER2)

#R PLOT
fig, ax1 = plt.subplots(1, 1, figsize=(5, 6))
ax1.errorbar(R_TmodP, R_diff, yerr=R_error, marker='+', label='$R$ Filtered', linestyle='none')
ax1.plot(time_t, F_t_R, color='r', linestyle='-', label=r'$F_R(t)$ Model')
ax1.set_xlabel('Time since initial observation (s)')
ax1.set_ylabel('$R_{\mathrm{V}1} - R_{\mathrm{ref}}$')
ax1.set_title(cluster_name + ' Variability of Star V1')
ax1.grid(True)
ax1.invert_yaxis()
ax1.legend(loc='lower right')
fig.show()

#V PLOT
fig, ax2 = plt.subplots(1, 1, figsize=(5, 6))
ax2.errorbar(V_TmodP, V_diff, yerr=V_error, marker='+',color='black', label='$V$ Filtered', linestyle='none')
ax2.plot(time_t, F_t_V, color='#32CD32', linestyle='-', label=r'$F_V(t)$ Model')
ax2.set_xlabel('Time since initial observation (s)')
ax2.set_ylabel('$V_{\mathrm{V}1} - V_{\mathrm{ref}}$')
ax2.set_title(cluster_name + ' Variability of Star V1')
ax2.grid(True)
ax2.invert_yaxis()
ax2.legend(loc='lower right')
plt.show()

#BOTH
fig, ax3 = plt.subplots(1, 1, figsize=(5, 6))
ax3.errorbar(R_TmodP, R_diff, yerr=R_error, marker='+', color='#1f77b4', ecolor='#1f77b4', label='$R$ Filtered', linestyle='none')
ax3.errorbar(V_TmodP, V_diff, yerr=V_error, marker='+', color='black', ecolor='black', label='$V$ Filtered', linestyle='none')
ax3.plot(time_t, F_t_V, color='#32CD32', linestyle='-', label=r'$F_V(t)$ Model')
ax3.plot(time_t, F_t_R, color='r', linestyle='-', label=r'$F_R(t)$ Model')
ax3.set_xlabel('Time since initial observation (s)')
ax3.set_ylabel('$M_{\mathrm{V1}} - M_{\mathrm{ref}}$')
ax3.set_title(cluster_name + ' Variability of Star V1')
ax3.grid(True)
ax3.invert_yaxis()
ax3.legend(loc='lower right')
plt.show()
