#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 10:15:15 2020

@author: niamhmurray
"""

import matplotlib.pyplot as plt
#import cartopy.crs as ccrs
import numpy as np
from netCDF4 import Dataset
import datetime
from datetime import date
from datetime import timedelta


def read_S2S(filename, fyear, iprint):

    '''
    Read in the ground variable data from the netCDF file.
    Input: name of file to read.
    Output: 
    :longitude  - degrees
    :latitude   - degrees
    :outdates   - calendar dates for data points
    :u          - zonal component of the wind (m/s)
    :v          - meridional component of the wind (m/s)
    :mslp       - mean sea level pressure (hPa)
    '''
    print()
    print('Reading file ',filename)
    data = Dataset(filename, 'r')
    if iprint == 1:
        print(data)
        print()
        print(data.dimensions)
        print()
        print(data.variables)
        print()
        
    ftime = data.variables['time'][:]
    alon = data.variables['longitude'][:]
    alat = data.variables['latitude'][:]
    u = data.variables['u10'][:,:,:,:]
    v = data.variables['v10'][:,:,:,:]
    mslp = data.variables['msl'][:,:,:,:]

    data.close()
    #
    # Time is in hours since 00UT, 1 Jan 1900.
    # Convert to timedelta format and then add to startcal to make a list of dates.
    # Note that you can add times in datetime and timedelta formats
    # which allows for leap years etc in time calculations.
    #
    startcal = datetime.datetime(1900, 1, 1)
    outdates = [startcal+timedelta(hours=ftimel) for ftimel in ftime]

    return alon, alat, outdates, u, v, mslp

def subset_field(alon, alat, lonpick, latpick):

    '''
    Find the indices of the grid point centred closest to chosen location.
    Input: 
    :alon       - longitude points
    :alat       - latitude points
    :lonpick    - longitude of chosen location
    :latpick    - latitude of chosen location
    Output:
    :intlon     - index of longitude for chosen point
    :intlat     = index of latitude for chosen point
    '''
    #
    # Using the fact that longitude and latitude are regularly spaced on grid.
    # Also, points ordered from north to south in latitude.
    # The indices (intlon, intlat) correspond to the centre of the closest grid-box 
    # to the chosen location.
    #
    dlon = alon[1]-alon[0]
    dlat = alat[1]-alat[0] # note that dlat is negative due to ordering of grid
    lonwest = alon[0]-0.5*dlon
    latnorth = alat[0]-0.5*dlat
    intlon = int(round((lonpick-lonwest)/dlon))
    intlat = int(round((latpick-latnorth)/dlat))
    print()
    print('Longitude of nearest grid box = ',alon[intlon])
    print('Latitude of nearest grid box = ',alat[intlat])
    
    return intlon, intlat

def plot_series(timarr, y, ylabel, mytitle):
    '''
    Plot the subset time series
    Inputs:
        timarr   - time array in datetime format
        y        - data time series
        ylabel   - string name for data
        mytitle  - plot title
    '''
    plt.figure()
    plt.plot(timarr,y,label=ylabel)
    plt.xlabel("Date")
    plt.ylim(0,17)
    plt.xticks(rotation=70)
    plt.ylabel(ylabel)
    plt.title(mytitle)
    plt.legend()
    plt.show()
    
def plot_histogram(array, bins, color, title, xlabel, ylabel):
    """plots a histogram for user choosen array
    Parameter: a: figure number
    Parameter: array: The array to be plotted
    Parameter: bins: The number of bins required for the plot 
    Parameter: color: the color required for the plot, input as string 
    Parameter: title: title of the plot, input as string 
    Parameter: xlabel: x-axis label, input as string 
    Parameter: ylabel: y-axis label, input as string 
    Parameter: xmin: The minimum value of x-axis
    Parameter: xmax: The maximum value of x-axis
    Returns: a Histogram plot plot """
    plt.hist(array,bins=bins,color=color, edgecolor='k')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim(0,100)
    plt.grid(linewidth="0.4")
    plt.axvline(np.percentile(array,50), color = 'k', linestyle = 'dashed', label = "50th Percentile")
    plt.axvline(np.percentile(array,75), color = 'r', linestyle = 'dashed', label = "75th Percentile")
    plt.axvline(np.percentile(array,98), color = 'g', linestyle = 'dashed', label = "98th Percentile")
    plt.legend()
    plt.show()
    
    median = np.median(array)
    upper_quartile = np.percentile(array, 75)
    v98 = np.percentile(array, 98)
    
    return median, upper_quartile, v98

def ssi_index(windarr, v98):
    """
    The function calculates and index for loss potential or Storm Serverity
    Index
    Where a loss is only assumed to have occured if the wind speeds exceed the 
    98th percentiale
    of wind speeds at the location. 
    The variable p is representative of the popultion density at the location 
    Parameter: windarr: an array containing a histroical data set of wind 
     speeds for the choosen locaiton 
    Parameter: v98: The 98th percentile of wind speeds at the choosen location 
    Returns: The potnetial loss over a given time for a choosen location 
    """
    index = np.zeros(len(windarr))
    
    for j in range(len(windarr)):
        if windarr[j] > v98:
            index[j] = (1*(windarr[j]/v98-1)**3)
        else:
            index[j] = 0
    return(index) 
    
def extract_series(fpath, fstem, lonpick, latpick, dstart, dend):
    '''
    High level function controlling extraction of runoff time series 
    for chosen location.
    Input: fpath, fstem determine the name of file to read
    :lonpick    - longitude of chosen location
    :latpick    - latitude of chosen location
    :dstart     - start date in datetime.date format
    :dend       - end date in datetime.date format
    Output: 
    :dayarr     - time in days since start
    :timarr     - time series in datetime format
    :windarr    - wind speed (m/s) time series at chosen location
    '''   
    #
    # Set end date and start date of required time series
    #
    dendp = dend+timedelta(days=1)
    tinterval = dendp-dstart
    ndays = tinterval.days
    #
    # Plot the data for the first date in the interval
    #
    fdate = dstart.strftime("%Y%m%d")
    fyear = dstart.year
    iprint = 0  # set to 1 to print variables on reading files; 0 for no print
    
    #Read the data
    filename = str(fpath+fstem+fdate+".nc")
    # Note that the str() function is included to ensure that these
    # variables are interpreted as character strings.
    
    alon, alat, outdates, u, v, mslp = read_S2S(filename, fyear, iprint)
    #
    # Find the indices of the grid box centred closest to the chosen location
    #
    intlon, intlat = subset_field(alon, alat, lonpick, latpick)
    
    # Setup arrays to save time series data
    #
    dayarr = np.arange(ndays)
    timarr = np.arange(np.datetime64(str(dstart)), np.datetime64(str(dendp)))
    windarr = np.zeros(ndays)
    #
    # Loop over months, reading files and saving daily data
    dcur = dstart
    for n in range(ndays):
        fdate = dcur.strftime("%Y%m%d")
        ilead = 5
        ensemble = 0
    #Read the data
        filename = str(fpath+fstem+fdate+'.nc')
    # Note that the str() function is included to ensure that these
    # variables are interpreted as character strings.
        alon, alat, outdates, u, v, mslp = read_S2S(filename, fyear, iprint)
    #
    # Save the data required from this time-point
        windarr[n] = np.sqrt(u[ilead, ensemble, intlat, intlon]**2 + v[ilead, ensemble, intlat, intlon]**2)
    #
    # Increment the date variable by one day
    #
        dcur=dcur+timedelta(days=1)


    return dayarr, timarr, windarr


if __name__ == '__main__':
    
    '''
    Main program script extracting time series from ERA5 data.
    '''
    #
    # Pick the location to extract the time series
# =============================================================================
#     # This coordinate is London.
#     lonpick = -0.1278
#     latpick = 51.5074
# =============================================================================

    #This coordinate is Paris
    lonpick = 2.3522
    latpick = 48.8566

# =============================================================================
#     #This coordinate is Hamburg 
#     lonpick = 9.9937
#     latpick = 53.5511
# =============================================================================

    
    # Select the start and end date required for the time series
    #
    dstart = datetime.date(1999, 12, 20)
    dend = datetime.date(1999, 12, 31)
    
    fpath = str('/Users/niamhmurray/Documents/Masters /Climate Services/Technical_Assingment/europe_NCEP_inst_energy/')
    fstem = str('eur_inst_energy_hc_')
    dayarr, timarr, windarr = extract_series(fpath, fstem, lonpick, latpick, dstart, dend)

    plot_series(timarr, windarr, 'wind Speed  (m/s)', 'Wind Speed Time Series \n for Paris from S2S Forecast during Storm Martin')
    
    #Use Plot Histogram to calculate v98 and show histogram
    median, upper_quartile, v98 = plot_histogram(windarr, 35, "palegreen", "Frequency of Windspeed \n in Paris 1999-2010 (DJF) from S2S Forecast \n with leadtime of 5 days","Wind Speed (m/s)", "Frequency")
    ssi_S2S_paris_martin = ssi_index(windarr, 8.09)
    acum_ssi_S2S_paris_martin = sum(ssi_S2S_paris_martin)
    