# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:55:44 2025

@author: lucie
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
import scipy.stats as stats
from scipy.stats import chi2
from matplotlib.patches import Ellipse

plt.ion()

folder_name = 'FINAL/M37'
cluster_name = 'M37'
name = 'M37'
file_path = '/Users/lucie/OneDrive/Documents/Uni Y4/Major Project/Calibration/NEW/'+folder_name+'/'  
 
excel_file = file_path+name+'_calibrated_magnitudes.xlsx'
F = pd.read_excel(excel_file, sheet_name=1)  #read V sheet of file
# Replace '--' with NaN and convert to float
F.replace('--', np.nan, inplace=True)  # Replaces '--' with NaN
F.dropna(subset=['PMRA', 'PMDEC', 'DEC'], inplace=True)  # Drop rows with missing values

PMRA = np.array(F['PMRA'], dtype=float)
PMDEC = np.array(F['PMDEC'], dtype=float)
index = np.array(F['Star Index'])
DEC = np.array(F['DEC'], dtype=float) 

fig, (ax1) = plt.subplots(1,1, figsize=(5,6))
ax1.scatter(PMRA, PMDEC, marker='+', label=' Cluster candidates with \n V image star indices \n for outliers')
ax1.set_xlabel(r'$\mu_\alpha$ [mas/yr]')
ax1.set_ylabel(r'$\mu_\delta$ [mas/yr]')
ax1.set_title(cluster_name+' Proper Motions')
ax1.grid(True)
ax1.legend()#location='upper right'
for i in range(len(PMRA)):
    x = PMRA[i]
    y = PMDEC[i]
    
    # Check if the star is in the specified region
    if not(-1<= x <=3 and -6.6<= y <=-4):
        ax1.text(x, y, str(index[i]), fontsize=8, color='black', ha='right')
        
        
#now we want to compare the membership probabilities of the cluster members
#using the mahalanobis distance (small = cluster member, large = field star)
def mahalanobis_distance(point, mean, inv_cov):
    diff = point-mean
    return np.sqrt(diff @ inv_cov @ diff.T) #@=matrix multiplication

PM = np.vstack([PMRA, PMDEC])
mean_pm = np.mean(PM, axis=1)
print('Mean PM: '+str(mean_pm))
cov_pm = np.cov(PM)
print("Covariance Matrix:\n", cov_pm)
cov_inv = np.linalg.inv(cov_pm)

#compute eigenvalues & eigenvectors (for ellipse orientation)
eigvals, eigvecs = np.linalg.eigh(cov_pm)
angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))  #rotation angle

#compute ellipse widths for 1σ and 2σ confidence regions
std_devs = np.sqrt(chi2.ppf([0.68, 0.95], df=2))
#print(std_devs)
widths = std_devs*np.sqrt(eigvals[0])
heights = std_devs*np.sqrt(eigvals[1])

#convert Mahalanobis distance to probability using chi2 CDF
distances = np.array([mahalanobis_distance(PM[:, i], mean_pm, cov_inv) for i in range(PM.shape[1])])
membership_probs = 1-stats.chi2.cdf(distances**2, df=2)
#print(membership_probs)
cluster_members = []
field_stars = []
for i,m in enumerate(membership_probs):
    if m >=0.8:
        cluster_members.append(i)
    elif m <0.8:
        field_stars.append(i)

PMDEC_cluster = []
PMRA_cluster = []
for idx in cluster_members:
    PMDEC_cluster.append(PMDEC[idx])
    PMRA_cluster.append(PMRA[idx])
    
PMDEC_field = []
PMRA_field = []
for idxs in field_stars:
    PMDEC_field.append(PMDEC[idxs])
    PMRA_field.append(PMRA[idxs])

fig1, (ax2) = plt.subplots(1,1, figsize=(5,6))
ax2.scatter(PMRA_cluster, PMDEC_cluster, color='r', marker='+', label='Cluster Members')
ax2.scatter(PMRA_field, PMDEC_field, color='b', marker='+', label='Field Stars')
ax2.set_xlabel(r'$\mu_\alpha \cos \delta$ [mas/yr]')
ax2.set_ylabel(r'$\mu_\delta$ [mas/yr]')
ax2.set_title(cluster_name+' Proper Motions with 1σ and 2σ Contours')
ax2.grid(True)
ax2.legend(loc='lower right')#location='upper right'
for i in range(len(PMRA)):
    x = PMRA[i]
    y = PMDEC[i]
    
    # Check if the star is in the specified region
    if not(0<= x <=4 and -6.5<= y <=-5):
        ax2.text(x, y, str(index[i]), fontsize=8, color='black', ha='right')
        
#add ellipses for 1σ and 2σ to plot
ellipse_1sigma = Ellipse(xy=mean_pm, width=widths[0], height=heights[0],
                         angle=angle, edgecolor='red', facecolor='none', linewidth=2, linestyle='solid')
ellipse_2sigma = Ellipse(xy=mean_pm, width=widths[1], height=heights[1],
                         angle=angle, edgecolor='blue', facecolor='none', linewidth=2, linestyle='dashed')

ax2.add_patch(ellipse_1sigma)
ax2.add_patch(ellipse_2sigma)

# Add labels to contours
ax2.text(mean_pm[0] + 1.1*widths[0], mean_pm[1], '1$\sigma$', color='red', fontsize=10, ha='center', va='center')
ax2.text(mean_pm[0] + 1.1*widths[1], mean_pm[1], '2$\sigma$', color='blue', fontsize=10, ha='center', va='center')
#add label of mean pm to plot
#btm left 0.05, 0.19
#top right 0.55, 0.85
ax2.text(0.55, 0.19, f'Mean PM: ({mean_pm[0]:.2f}, {mean_pm[1]:.2f})', transform=ax2.transAxes, fontsize=10, ha='left', va='top', 
         bbox=dict(facecolor='white', edgecolor='grey', boxstyle='round,pad=0.5'))