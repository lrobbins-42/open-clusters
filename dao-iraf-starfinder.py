# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 12:25:29 2024

@author: lucie
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:36:58 2024

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
    plt.imshow(B2D_value, origin='lower', cmap='inferno', aspect='equal')
    ax.invert_yaxis()
    plt.title(cluster_name + ' Background', y=1.05)
    cbar = plt.colorbar(ax=ax, orientation='vertical', pad=0.08, fraction=0.046)
    cbar.set_label('Intensity', fontsize=10) 
    #plt.colorbar()
    plt.xlabel('RA')
    plt.ylabel('DEC')
    
def StarIdentifierDao(cluster_name):
    '''Circles pixels centres identified as stars and plots using DAOStarFinder''' 
    
    positions = np.array(list(zip(dao['xcentroid'], dao['ycentroid'])))
    #print("DAO Star Positions (pixels):", positions)
    
    # Convert pixel coordinates to world coordinates (RA, Dec)
    DAOworld_coords = list(map(tuple, wcs_cluster.all_pix2world(positions, 0)))
    #print("DAO World Coordinates (RA, Dec):", DAOworld_coords)
    DAO_RA, DAO_DEC = zip(*DAOworld_coords)
    magnitudes = np.array(list(zip(dao['mag'])))
    #print(magnitudes)
    
    dao_dict = {'RA': DAO_RA, 'DEC': DAO_DEC, FILTER: magnitudes.flatten()}
    dao_df = pd.DataFrame(dao_dict)

    excel_file = 'NGC7686_dao_magnitudes_'+FILTER+'.xlsx'
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        dao_df.to_excel(writer, index=False, sheet_name='Sheet 1', float_format='%.4f')

        worksheet1 = writer.sheets['Sheet 1']
      
        for col_num, value in enumerate(dao_df.columns.values):
            worksheet1.column_dimensions[chr(65 + col_num)].width = 12
        print(f'\nSummary data has been written to {excel_file}')
    
    # Circles Star Locations
    apertures = CircularAperture(positions, r=3)
    
    # Normalize and plot the image
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection=wcs_cluster)
    plt.grid()
    #overlay = ax.get_coords_overlay('icrs')
    #overlay.grid(color='black', ls='dotted')
    #norm = ImageNormalize(stretch=SqrtStretch())
    plt.title(cluster_name+" DAO Star Identifier: Sigma = "+str(sigma_T)+"" , loc='left')
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

def StarIdentifierIRAF(cluster_name):
    '''Circles pixels centres identified as stars and plots using IRAFStarFinder'''    
    positions = np.array(list(zip(IRAF['xcentroid'], IRAF['ycentroid'])))
    #print("IRAF Star Positions (pixels):", positions)
    
    # Convert pixel coordinates to world coordinates (RA, Dec)
    IRAFworld_coords = list(map(tuple, wcs_cluster.all_pix2world(positions, 0)))
    #print("IRAF World Coordinates (RA, Dec):", IRAFworld_coords)
    
    # Convert pixel coordinates to world coordinates (RA, Dec)
    IRAFworld_coords = list(map(tuple, wcs_cluster.all_pix2world(positions, 0)))
    #print("DAO World Coordinates (RA, Dec):", DAOworld_coords)
    IRAF_RA, IRAF_DEC = zip(*IRAFworld_coords)
    magnitudes = np.array(list(zip(IRAF['mag'])))
    #print(magnitudes)
    
    IRAF_dict = {'RA': IRAF_RA, 'DEC': IRAF_DEC, FILTER: magnitudes.flatten()}
    IRAF_df = pd.DataFrame(IRAF_dict)

    excel_file = 'NGC7686_IRAF_magnitudes_'+FILTER+'.xlsx'
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        IRAF_df.to_excel(writer, index=False, sheet_name='Sheet 1', float_format='%.4f')

        worksheet1 = writer.sheets['Sheet 1']
      
        for col_num, value in enumerate(IRAF_df.columns.values):
            worksheet1.column_dimensions[chr(65 + col_num)].width = 12
        print(f'\nSummary data has been written to {excel_file}')
    
    # Circles Star Locations
    apertures = CircularAperture(positions, r=3)
    
    # Normalize and plot the image
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection=wcs_cluster)
    plt.grid()
    #overlay = ax.get_coords_overlay('icrs')
    #overlay.grid(color='black', ls='dotted')
    #norm = ImageNormalize(stretch=SqrtStretch())
    plt.title(cluster_name+" IRAF Star Identifier: Sigma = "+str(sigma_T)+"" , loc='left')
    mid = np.median(image)
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.imshow(image, cmap='Greys', origin='lower', vmin=mid-(0.3*mid), vmax=mid+(0.5*mid))
    ax.invert_yaxis()
    cbar = plt.colorbar(ax=ax, orientation='vertical', pad=0.05, fraction=0.046)
    cbar.set_label('Intensity', fontsize=10) 
    #plt.colorbar()
    
    # Plot the apertures on the image
    apertures.plot(color='red', lw=1.5, alpha=1)
    
     # Add labels for each detected star
    for i, (x, y) in enumerate(positions):
        # Add a label at the (x, y) coordinates
        ax.text(x + 5, y + 5, str(i+1), color='white', bbox=dict(facecolor='black', alpha=0.2, pad=0.2), fontsize=10, fontweight='bold')
        
    return IRAFworld_coords
        
def Both(cluster_name):
    '''Plots DAOPHOT and IRAF star candidates on one plot'''
    dao_positions = np.array(list(zip(dao['xcentroid'], dao['ycentroid'])))
    IRAF_positions = np.array(list(zip(IRAF['xcentroid'], IRAF['ycentroid'])))
    # Circles Star Locations
    dao_apertures = CircularAperture(dao_positions, r=7)
    IRAF_apertures = CircularAperture(IRAF_positions, r=3)
    
    # Normalize and plot the image
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot(projection=wcs_cluster)
    plt.grid()
    #overlay = ax.get_coords_overlay('icrs')
    #overlay.grid(color='black', ls='dotted')
    #norm = ImageNormalize(stretch=SqrtStretch())
    plt.title(cluster_name+" DAO(blue) and IRAF(red) Star Identifiers: Sigma = "+str(sigma_T)+"" , loc='left')
    mid = np.median(image)
    plt.imshow(image, cmap='Greys', origin='lower', vmin=mid-(0.3*mid), vmax=mid+(0.5*mid))
    ax.invert_yaxis()
    plt.xlabel('RA')
    plt.ylabel('DEC')
    cbar = plt.colorbar(ax=ax, orientation='vertical', pad=0.05, fraction=0.046)
    cbar.set_label('Intensity', fontsize=10) 
    #plt.colorbar()
    
    # Plot the apertures on the image
    dao_apertures.plot(color='blue', lw=1.5, alpha=1)
    IRAF_apertures.plot(color='red', lw=1.5, alpha=1)
    
    for i, (x, y) in enumerate(dao_positions):
        # Add a label at the (x, y) coordinates
        ax.text(x + 5, y + 5, str(i+1), color='white', bbox=dict(facecolor='black', alpha=0.2, pad=0.2), fontsize=10, fontweight='bold')
    
     
    
def DATABASE(RA,DEC):        
    ''' Finds SIMBAD Coordinates of Target Cluster and prints them in a table'''     
    customSimbad = Simbad()
    customSimbad.TIMEOUT = 120  # Adjust timeout in case of large queries
    #customSimbad.remove_votable_fields('coordinates')  # Remove unnecessary fields
    customSimbad.add_votable_fields('flux(B)', 'flux(V)')  # Add B and V magnitudes

    c = SkyCoord(ra=RA, dec=DEC, frame='icrs', unit=(u.hourangle, u.deg))
    
    r = 2 * u.arcminute #searches database for stars within 'r' of coordinates inserted into function
    ra_deg=c.ra.deg
    dec_deg=c.dec.deg
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
    excel_file = 'NGC7686_'+FILTER+'_simbad_magnitudes.xlsx'
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        df.to_excel(writer, index=False, sheet_name='Sheet 1', float_format='%.4f')

        worksheet1 = writer.sheets['Sheet 1']
      
        for col_num, value in enumerate(df.columns.values):
            worksheet1.column_dimensions[chr(65 + col_num)].width = 12
        print(f'\nSummary data has been written to {excel_file}')


plt.ion()

#change name, filter and exp below when changing target image
#DAO/IRAF will return magnitudes based on the filter selected

#TARGET IMAGE: NGC7686-0001_G_30.fits
name = 'NGC 7686'
FILTER = 'V'
exp = '30'
version = '0001'
#image_file = 'm34_30s_'+FILTER+'_1.fit'
image_file = 'NGC7686-'+version+'_'+FILTER+'_'+exp+'.fits'
#format figure titles with cluster_name
cluster_name = name+' '+exp+'/'+FILTER

header_data_unit_list = fits.open(image_file)

image = header_data_unit_list[0].data
header = header_data_unit_list[0].header

#print(header)

wcs_cluster = WCS(header)
print(wcs_cluster)


'''Used for defining a threshold value for both algorithms
above which sources are detected'''
sigma_T = float(input('sigma value for identifying stars: ='))  

''' global mean. median and standard deviation for image'''
mean, median, std = sigma_clipped_stats(image, sigma=7)
#print((mean, median, std))

'''Background calculation using Background2D'''
B2D = Background2D(image,101)
B2D_value=B2D.background
#print(B2D_value)


'''Runs DAOPHOT algorithm and subtracts Background2D'''
daofind = DAOStarFinder(fwhm=5.0, threshold=sigma_T*std, brightest=(15))#, min_separation=3,)
dao = daofind(image - B2D_value)
for col in dao.colnames:
    dao[col].info.format = '%.8g' 
#print(dao)



'''Runs IRAF algorithm and subtracts Background2D'''
IRAFfind = IRAFStarFinder(fwhm=5.0,threshold=sigma_T*std, brightest=(15))#, , min_separation=3)
IRAF = IRAFfind(image-B2D_value)
for col in IRAF.colnames:
    IRAF[col].info.format = '%.8g' 
#print(IRAF)


Image(cluster_name)
Background(cluster_name)

DAO_cb = StarIdentifierDao(cluster_name)
IRAF_cb = StarIdentifierIRAF(cluster_name)

DAO_RA, DAO_DEC = zip(*DAO_cb)
IRAF_RA, IRAF_DEC = zip(*IRAF_cb)
Both(cluster_name)

#DATABASE('02h42m07.4s' ,'+42d43m19s') #M34
DATABASE('23h29m41.6s' ,'+49d10m29s') #NGC 7686





