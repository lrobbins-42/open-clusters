# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 16:20:44 2025

@author: lucie
"""
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
import photutils
import astroquery
from photutils.detection import DAOStarFinder
from photutils.detection import IRAFStarFinder
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
from astropy import coordinates
import astropy.units as u
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles

from photutils.background import Background2D, MedianBackground, MMMBackground
import numpy as np
import pandas as pd
from astropy.coordinates import match_coordinates_sky
from photutils.aperture import CircularAperture, aperture_photometry

'''wcs_input_dict = {
    'CTYPE1': 'RA---TAN', 
    'CUNIT1': 'deg', 
    'CDELT1': -0.0002777777778, 
    'CRPIX1': 1, 
    'CRVAL1': 337.5202808, 
    'NAXIS1': 1024,
    'CTYPE2': 'DEC--TAN', 
    'CUNIT2': 'deg', 
    'CDELT2': 0.0002777777778, 
    'CRPIX2': 1, 
    'CRVAL2': -20.833333059999998, 
    'NAXIS2': 1024
}
wcs_cluster_dict = WCS(wcs_input_dict)'''
def Image(cluster_name):
    '''Plots original image'''
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection=wcs_cluster)
    mid = np.median(image)
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='black', ls='dotted')
    overlay[0].set_axislabel(' ')
    overlay[1].set_axislabel(' ')
    plt.imshow(image, origin='lower', cmap='Greys', aspect='equal', vmin=mid-(0.3*mid), vmax=mid+(0.5*mid))
    ax.invert_yaxis()
    plt.title(cluster_name+' Starfield Image', y=1.05)
    cbar = plt.colorbar(ax=ax, orientation='vertical', pad=0.08, fraction=0.046)
    cbar.set_label('Intensity', fontsize=10) 
    #plt.colorbar()
    plt.xlabel('RA')
    plt.ylabel('DEC')
    
def Background(cluster_name):
    '''Plots the background image'''
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection=wcs_cluster)
    overlay = ax.get_coords_overlay('icrs')
    overlay.grid(color='white', ls='dotted')
    overlay[0].set_axislabel(' ')
    overlay[1].set_axislabel(' ')
    plt.imshow(background_image, origin='lower', cmap='inferno', aspect='equal')
    ax.invert_yaxis()
    plt.title(cluster_name + ' Background', y=1.05)
    cbar = plt.colorbar(ax=ax, orientation='vertical', pad=0.08, fraction=0.046)
    cbar.set_label('Intensity', fontsize=10) 
    #plt.colorbar()
    plt.xlabel('RA')
    plt.ylabel('DEC')
    
def perform_aperture_photometry(image, positions, fwhm, exp):
    aperture_radius = 3*fwhm
    apertures = CircularAperture(positions, r=aperture_radius)
    phot_table = aperture_photometry(image, apertures)
    #GAIN = 2.3900001049041748
    GAIN = 2.39
    aperture_sum = phot_table['aperture_sum']
    N_counts = (GAIN*aperture_sum)/exp
    
    return N_counts
    
def StarIdentifierDao(name, cluster_name, exp, dao, image, image_bg_subtracted, wcs_cluster, fwhm):
    '''Circles pixels centres identified as stars and plots using DAOStarFinder''' 
    exp = exp
    positions = np.array(list(zip(dao['xcentroid'], dao['ycentroid'])))
    #print("DAO Star Positions (pixels):", positions)
    
    # Convert pixel coordinates to world coordinates (RA, Dec)
    DAOworld_coords = list(map(tuple, wcs_cluster.all_pix2world(positions, 0)))
    #print("DAO World Coordinates (RA, Dec):", DAOworld_coords)
    DAO_RA, DAO_DEC = zip(*DAOworld_coords)
    
    N_counts = perform_aperture_photometry(image_bg_subtracted, positions, fwhm, exp)
    #extract counts (N)
    #compute magnitudes and errors using Poisson stats
    inst_mags = -2.5*np.log10(N_counts)
    sigma_N = np.sqrt(N_counts)
    inst_mag_errors = (2.5/np.log(10))*(sigma_N/N_counts)
    
    phot_df = pd.DataFrame({'RA': DAO_RA,'DEC': DAO_DEC,FILTER: inst_mags,FILTER+' Error': inst_mag_errors})

    excel_file = name+'_dao_magnitudes_'+FILTER+'.xlsx'
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        phot_df.to_excel(writer, index=False, sheet_name='Sheet 1', float_format='%.4f')

        worksheet1 = writer.sheets['Sheet 1']
      
        for col_num, value in enumerate(phot_df.columns.values):
            worksheet1.column_dimensions[chr(65 + col_num)].width = 12
        print(f'\nSummary data has been written to {excel_file}')
    
    # Circles Star Locations
    apertures = CircularAperture(positions, r=3)
    
    # Normalize and plot the image
    plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection=wcs_cluster)
    plt.grid()
    #overlay = ax.get_coords_overlay('icrs')
    #overlay.grid(color='black', ls='dotted')
    #norm = ImageNormalize(stretch=SqrtStretch())
    plt.title(cluster_name+" DAO Star Identifier: Sigma = "+f'{sigma_T:.2f}'+"" , loc='left')
    mid = np.median(image)
    plt.imshow(image, cmap='Greys', origin='lower', vmin=mid-(0.3*mid), vmax=mid+(0.5*mid))
    ax.invert_yaxis()
    plt.xlabel('RA')
    plt.ylabel('DEC')
    cbar = plt.colorbar(ax=ax, orientation='vertical', pad=0.05, fraction=0.046)
    cbar.set_label('Intensity', fontsize=10) 
    #plt.colorbar()
    
    # Plot the apertures on the image
    apertures.plot(color='blue', lw=1.5, alpha=1)
    
      # Add labels for each detected star
    for i, (x, y) in enumerate(positions):
        # Add a label at the (x, y) coordinates
        ax.text(x + 5, y + 5, str(i+1), color='white', bbox=dict(facecolor='black', alpha=0.2, pad=0.2), fontsize=10, fontweight='bold')

    return DAOworld_coords

    
def DATABASE(name,RA,DEC,search_width):        
    ''' Finds SIMBAD Coordinates of Target Cluster and prints them in a table'''     
    customSimbad = Simbad()
    customSimbad.TIMEOUT = 120  # Adjust timeout in case of large queries
    #customSimbad.remove_votable_fields('coordinates')  # Remove unnecessary fields
    customSimbad.add_votable_fields('flux(B)', 'flux(V)')  # Add B and V magnitudes

    c = SkyCoord(ra=RA, dec=DEC, frame='icrs', unit=(u.hourangle, u.deg))
    
    r = search_width*u.arcminute #searches database for stars within 'r' of coordinates inserted into function
    #ra_deg=c.ra.deg
    #dec_deg=c.dec.deg
    result_table = customSimbad.query_region(c, radius=r)
    print(result_table)
    
    valid_data = []

    for row in result_table:
        # Check if the 'B' and 'V' magnitudes are valid (not '--')
        B_mag = row['FLUX_B']
        V_mag = row['FLUX_V']
        
        if B_mag != '--' and V_mag != '--':  
            # Convert RA and DEC to degrees (if needed, depending on the format)
            ra_result_deg = SkyCoord(row['RA'], row['DEC'], unit=(u.hourangle, u.deg)).ra.deg
            dec_result_deg = SkyCoord(row['RA'], row['DEC'], unit=(u.hourangle, u.deg)).dec.deg
            
            # Append the valid data to the list
            valid_data.append({
                'ID': row['MAIN_ID'],
                'RA (deg)': ra_result_deg,
                'DEC (deg)': dec_result_deg,
                'B': B_mag,
                'V': V_mag
            })


    df = pd.DataFrame(valid_data)
    excel_file = name+FILTER+'_simbad_magnitudes.xlsx'
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Sheet 1', float_format='%.4f')

        worksheet1 = writer.sheets['Sheet 1']
      
        for col_num, value in enumerate(df.columns.values):
            worksheet1.column_dimensions[chr(65 + col_num)].width = 12
        print(f'\nSummary data has been written to {excel_file}')

    return df

plt.ion()

#change name, filter and exp below when changing target image
#DAO will return magnitudes based on the filter selected

#TARGET IMAGE: NGC7686-0001_G_30.fits
#m34_30s_b_1.fit
name = 'NGC744'
FILTER = 'B'
exp = '60s'
EXPOSURE_TIME = 60
version = '0004'
#image_file = name+'-'+version+'_'+FILTER+'_'+exp+'.fits'
image_file = name+'-'+version+'_'+exp+'_'+FILTER+'.fits'
#image_file = 'm34_30s_v_1.fit'
#format figure titles with cluster_name
cluster_name = name+' '+exp+'/'+FILTER

header_data_unit_list = fits.open(image_file)

image = header_data_unit_list[0].data
header = header_data_unit_list[0].header

#print(header)

wcs_cluster = WCS(header)
print(wcs_cluster)


#sigma_T = float(input('sigma value for identifying stars: ='))  
#background_mean = np.median(background_image)
#background_std = np.std(background_image)
#B2D_value=B2D.background
#print(B2D_value)


'''Background calculation using Background2D RMS noise = 
variation in pixel values caused by random noise sources in the background'''
B2D = Background2D(image,151)
background_image = B2D.background
#mean_rms = np.mean(B2D.background_rms)
#std_rms = np.std(B2D.background_rms)
#threshold = mean_rms + (5*std_rms) #5 sigma above the background fluctuation
#background_mean = np.median(background_image)
#background_std = np.std(background_image)
#threshold = background_mean + (5*background_std) 
threshold = 100
print("threshold: ",threshold)  


#FWHM=10
fwhm=12
'''Runs DAOPHOT algorithm and subtracts Background2D'''
#SPECIFY A SIGMA_RADIUS OF 3 RATHER THAN THE PRESUMED VALUE OF 1.5
daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold, brightest=(35))#, min_separation=3,)
image_bg_subtracted = image - background_image
dao = daofind(image_bg_subtracted)
sigma_T = fwhm/2.355
#for col in dao.colnames:
 #   dao[col].info.format = '%.8g' 
#print(dao)

Image(cluster_name)
Background(cluster_name)

DAO_cb = StarIdentifierDao(name, cluster_name, EXPOSURE_TIME, dao, image, image_bg_subtracted, wcs_cluster, fwhm)
DAO_RA, DAO_DEC = zip(*DAO_cb)
DAO_RA = np.array(DAO_RA)
DAO_DEC = np.array(DAO_DEC)
#print(len(DAO_RA))
#print(len(DAO_DEC))

search_width = 7 #CHANGE BASED ON SIMBAD IMAGE 
#df = DATABASE(name,'02h42m07.4s' ,'+42d43m19s', search_width) #M34
#df = DATABASE(name, '23h29m41.3s', '+49d10m12s', search_width) #NGC 7686
df = DATABASE(name, '01h58m36.5s', '+55d28m23s', search_width) #NGC 744
#df = DATABASE(name, '05h36m20.2s', '+34d08m06s', search_width) #M36
SIM_RA = df['RA (deg)'].values
#print(SIM_RA)
SIM_DEC = df['DEC (deg)'].values

mid = np.median(image)
plt.figure(figsize=(10, 10))
ax = plt.subplot(projection=wcs_cluster)
plt.imshow(image, origin='lower', cmap='Greys', vmin=mid-(0.3*mid), vmax=mid+(0.5*mid))
ax.scatter(DAO_RA, DAO_DEC, transform=ax.get_transform('world'), s=50, edgecolor='blue', facecolor='none', label='DAOPHOT')
ax.scatter(SIM_RA, SIM_DEC, transform=ax.get_transform('world'), s=50, edgecolor='red', facecolor='none', label='SIMBAD')
ax.invert_yaxis()
plt.legend()
plt.show()

#from astropy.io import fits
#hdul = fits.open(image_file)
#hdr = hdul[0].header
#print(hdr)