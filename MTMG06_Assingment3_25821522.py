#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 09:49:38 2020

@author: niamhmurray
"""

import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
from netCDF4 import num2date, date2num, Dataset as NetCDFFile
from mpl_toolkits.basemap import Basemap 
import scipy.stats as stats


def calculate_area(lats):
    """
    Returns the area within a grid box dependant on the latitude of the gird box
    Parameter: lats: reguired latitudinal areas
    output: Area depedant on Latitude
    """

    area = np.sqrt(np.cos(lats*np.pi/180))
    return area

def data_covariance_matrix():
    """
    Converts a 3D raw data set into a 2D data set then sclaes it depending on 
    grid box area
    Returns: x: 2D unsclaed data matrix
    Returns: x_tilda: 2D scaled Data Matirx 
    Returns: s: Covariance Matrix of the sclaed data matrix
    Returns: f: Diagonal scaling matrix (Latitudinal dependance) 
    Returns: m: Climatological Mean as a list
    """
    columns = len(lons) * len(lats)
    rows = len(a_p[:,0,0])
    x = np.zeros((rows, columns))
    x_tilda = np.zeros_like(x)
    f = np.zeros((columns,columns))
    m = np.zeros(columns)
    area = calculate_area(lats)
    for ii in range (rows):
        for jj in range(len(lats)):
            for kk in range(len(lons)):
                f[(jj*len(lons))+kk, (jj*len(lons))+kk] = area[jj]
                if m[jj*len(lons)+kk] == 0:
                    m[jj*len(lons)+kk] = mu[jj,kk]
                x[ii,(jj*len(lons))+kk] = a_p[ii,jj,kk]
        x_tilda[ii] = np.dot(x[ii]- m, f)
    
    s = np.cov(x_tilda, y = None, rowvar = False)
    
    return x, x_tilda, s, f, m, columns

def Contour_Plot(x, title):
    """
    Plots pressure countours for the north of Europe
    Parameter: x: the pressure Array to be plotted
    Parameter: title: (Input as a string), Desired title for the plot 
    Returns: Pressure Contour Basemap plot
    """
    fig_mu = plt.figure()
    ax_mu = fig_mu.add_subplot(111)
    ax_mu.set_title(title)
    m_mu = Basemap(width = 9000000, height = 5500000, resolution = 'l', projection = 'laea', lat_ts = 50, lat_0 = 50, lon_0 = -10)
    m_mu.drawcoastlines()
    m_mu.drawparallels(np.arange(25, 75, 10))
    m_mu.drawmeridians(np.arange(-60, 60, 10))
    grid = np.meshgrid(lons, lats)
    c,d = m_mu(*grid)
    clevs = np.arange(1002, 1024, 3)
    cs = m_mu.contourf(c,d, x, clevs, extend = 'both', cmap = 'RdYlBu_r' )
    cbar = m_mu.colorbar(cs, location = 'bottom')
    fig_mu.show()
    
def Plot_Depression(x, title):
    """
    Plots pressure countours for the north of Birthish isles
    Parameter: x: the pressure Array to be plotted
    Parameter: title: (Input as a string), Desired title for the plot 
    Returns: Pressure Contour Basemap plot of british isles
    """
    fig_mu = plt.figure()
    ax_mu = fig_mu.add_subplot(111)
    ax_mu.set_title(title)
    m_mu = Basemap(llcrnrlon=-20.5,llcrnrlat=45,urcrnrlon=15.,urcrnrlat=57., resolution = 'l', projection = 'laea', lat_0 = 45, lon_0 = -20)
    m_mu.drawcoastlines()
    m_mu.drawparallels(np.arange(45, 60, 4))
    m_mu.drawmeridians(np.arange(-20, 15, 4))
    grid = np.meshgrid(lons, lats)
    c,d = m_mu(*grid)
    clevs = np.arange(1014, 1024, 1)
    cs = m_mu.contourf(c,d, x, clevs, extend = 'both', cmap = 'RdYlBu_r' )
    cbar = m_mu.colorbar(cs, location = 'bottom')
    fig_mu.show()
    
def PCA_Contour_plot(p, title):
    """
    Plots PCA countours for the north of Europe
    Parameter: x: the PCA Loading Array to be plotted
    Parameter: title: (Input as a string), Desired title for the plot 
    Returns: PCA Contour Basemap plot
    """
    fig_mu = plt.figure()
    ax_mu = fig_mu.add_subplot(111)
    ax_mu.set_title(title)
    m_mu = Basemap(width = 9000000, height = 5500000, resolution = 'l', projection = 'laea', lat_ts = 50, lat_0 = 50, lon_0 = -10)
    m_mu.drawcoastlines()
    m_mu.drawparallels(np.arange(25, 75, 10))
    m_mu.drawmeridians(np.arange(-60, 60, 10))
    grid = np.meshgrid(lons, lats)
    c,d = m_mu(*grid)
    clevs = np.arange(-0.1, 0.1, 0.001)
    cs = m_mu.contourf(c,d, p, clevs, extend = 'both', cmap = 'RdYlBu_r' )
    cbar = m_mu.colorbar(cs, location = 'bottom')
    fig_mu.show()
                   
def PCA_Contour_plot_scaled(p, title):
    """
    Plots PCA countours scaled to pressure anaomaly for the north of Europe
    Parameter: x: the Pressure Array to be plotted
    Parameter: title: (Input as a string), Desired title for the plot 
    Returns: Pressure anomaly Contour Basemap plot
    """
    fig_mu = plt.figure()
    ax_mu = fig_mu.add_subplot(111)
    ax_mu.set_title(title)
    m_mu = Basemap(width = 9000000, height = 5500000, resolution = 'l', projection = 'laea', lat_ts = 50, lat_0 = 50, lon_0 = -10)
    m_mu.drawcoastlines()
    m_mu.drawparallels(np.arange(25, 75, 10))
    m_mu.drawmeridians(np.arange(-60, 60, 10))
    grid = np.meshgrid(lons, lats)
    c,d = m_mu(*grid)
    clevs = np.arange(-3, 3, 0.05)
    cs = m_mu.contourf(c,d, p, clevs, extend = 'both', cmap = 'RdYlBu_r' )
    cbar = m_mu.colorbar(cs, location = 'bottom')
    fig_mu.show()

def plot_time_series (array1, array2, c, title, xlabel, ylabel):
    """plots a time series for user choosen arrays
    Parameter: a: figure number
    Parameter: array1: The x values required for the time series
    Parameter: array2: The y values required for time series
    Parameter: c: the color required for the timeseries, input as string 
    Parameter: title: title of the plot, input as string 
    Parameter: xlabel: x-axis label, input as string 
    Parameter: ylabel: y-axis label, input as string 
    Parameter: ymin: The minimum value of y-axis
    Parameter: ymax: The maximum value of y-axis
    Returns: a time series plot 
    """
    plt.figure()
    plt.plot(array1, array2, color = c)
    plt.title(title)
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
 
def truncated_Depression_plot(trunc, title):
    """
    """
    coefficents = u[6653,:]
    e_scaled_coeff = np.zeros((931,trunc))
    for i in range(0,trunc):
        e_scaled_coeff[:,[i]] =(coefficents[i]*scaled_eigen_vectors[:,[i]])

    sum_scale = np.zeros((931))
    for g in range(0,931):
        sum_scale[g] = sum(e_scaled_coeff[g])

    m_scale = (m + sum_scale)
    trunceigen =np.resize(m_scale, [19,49])

    Plot_Depression(trunceigen, title)
    
# =============================================================================
# def compostite_loop_pos(variable):
#     """
#     """
#     array_all_positive = np.empty((36,76))
#     for j in range ((76)):
#         for i in range ((36)):
#             if re_snao[i,j] > 0 :
#                 array_all_positive[i,j] = variable[i,j]
#             else:
#                 array_all_positive[i,j] = ('nan')
#             
#     return array_all_positive 
# 
# def compostite_loop_neg(variable):
#     """
#     """
#     array_all_negative = np.empty((36,76))
#     for j in range ((76)):
#         for i in range ((36)):
#             if re_snao[i,j] < 0 :
#                 array_all_negative[i,j] = variable[i,j]
#             else:
#                 array_all_negative[i,j] = ('nan')
#     
#     return array_all_negative 
# =============================================================================

def composition_plot(title, array, low, upper, interval):
    """
    """
    fig_comp = plt.figure()
    ax_comp = fig_comp.add_subplot(111)
    ax_comp.set_title(title)
    m_comp = Basemap(width = 9000000, height = 5500000, resolution = 'l', projection = 'laea', lat_ts = 50, lat_0 = 50, lon_0 = -10)
    m_comp.drawcoastlines()
    m_comp.drawparallels(np.arange(25,75,10))
    m_comp.drawmeridians(np.arange(-60,60,10))
    grid = np.meshgrid(p_lons, p_lats)
    x,y = m_comp(*grid)
    clevs = np.arange(low, upper, interval)
    cs = m_comp.contourf(x,y, array, clevs, extend = 'both', cmap = 'RdYlBu_r')
    cbar = m_comp.colorbar(cs, location = 'bottom')
    fig_comp.show()

def composite_analysis(variable):
    """
    """
    SNAO = u[:,0]
    #The overlap of Pressure and Variable data
    snao_array =SNAO[3100:6882]
    data = np.resize(variable, (3968, 36*76))
    #Remove the extra years from data so they line up 
    data_2 = data[0:3782, 0:2736]
    for i in range(3782):
        for j in range(2736):
            if data_2[i,j]<-100:
                data[i,j]= np.NAN
    positive = snao_array>0
    negative = snao_array<0
    #pull out overlap of data
    positive_snao = data_2[positive]
    negative_snao = data_2[negative]
    #find mean when SNAO is positive an reshape
    positive_mean = np.mean(positive_snao, axis = 0)
    positve_mean_plot = np.resize(positive_mean, (36, 76))
    #find mean when SNAO is negative an reshape
    negative_mean = np.mean(negative_snao, axis = 0)
    negative_mean_plot = np.resize(negative_mean, (36, 76))
    #Find the Difference between positive SNAO and negative SNAO means
    mean_diff = positve_mean_plot - negative_mean_plot
    
    return mean_diff, positve_mean_plot, negative_mean_plot, positive_mean, negative_mean
    
    
def stipple_composition_plot(variable, stippling_variable, title):
    """
    """
    fig_comp = plt.figure()
    ax_comp = fig_comp.add_subplot(111)
    ax_comp.set_title(title)
    m_comp = Basemap(width = 9000000, height = 5500000, resolution = 'l', projection = 'laea', lat_ts = 50, lat_0 = 50, lon_0 = -10)
    m_comp.drawcoastlines()
    m_comp.drawparallels(np.arange(25,75,10))
    m_comp.drawmeridians(np.arange(-60,60,10))
    grid = np.meshgrid(p_lons, p_lats)
    x,y = m_comp(*grid)
    dp_clevs = [-2, -1, -0.5, 0, 0.5, 1, 2]
    cs_cd = m_comp.contourf(x,y, variable, dp_clevs, extend = 'both', cmap = 'BrBG')
    cbar_cd = m_comp.colorbar(cs_cd, location = 'right')
    pv_levs = [0, 0.05, 1]
    m_comp.contourf(x,y, stippling_variable, pv_levs, extend = 'neither', hatches = ['.', None])
    fig_comp.show()
        
    
# open NetCDF file.
input_file = str("/Users/niamhmurray/Documents/Masters /Stats /Assingments/Assingment 3/MTMG06_ass3_slp.nc")
nc = NetCDFFile(input_file)



#Exploratory Analysis of the data
print(nc)
print(nc.dimensions)
print(nc.variables)

ncvar = 'msl'

#Lon and Lat Coordinates extracted from the NetCDF File
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

#Create a pressure array
a_p = nc.variables[ncvar][:]

#Read in time data and convert into an approprite unit
raw_times = nc.variables['time'][:]
times = num2date(raw_times, units= nc.variables['time'].units, calendar = nc.variables['time'].calendar)

#Close the NetCDF File 
nc.close()

#Read in Precipitation data
p_ffn = str("/Users/niamhmurray/Documents/Masters /Stats /Assingments/Assingment 3/MTMG06_ass3_precip.nc")
p_nc = NetCDFFile(p_ffn)

# =============================================================================
# #Exploratory Analysis of data 
# print(p_nc)
# print(p_nc.dimensions)
# print(p_nc.variables)
# =============================================================================

#Lon and Lat Coordinates extracted from the NetCDF File
p_lons = p_nc.variables['lon'][:]; n_p_lon = len(p_lons)
p_lats = p_nc.variables['lat'][:]; n_p_lat = len(p_lats)

#Read in time date and convert to a usable structure 
a_pr = p_nc.variables['precip'][:]
p_rtimes = p_nc.variables['time'][:]
p_times = num2date(p_rtimes, units = p_nc.variables['time'].units , calendar = p_nc.variables['time'].calendar)

#Close File 
p_nc.close()

#Read in Precipitation data
t_ffn = str("/Users/niamhmurray/Documents/Masters /Stats /Assingments/Assingment 3/MTMG06_ass3_temp.nc")
t_nc = NetCDFFile(t_ffn)

# =============================================================================
# #Exploratory Analysis of data 
# print(t_nc)
# print(t_nc.dimensions)
# print(t_nc.variables)
# =============================================================================

#Lon and Lat Coordinates extracted from the NetCDF File
t_lons = t_nc.variables['lon'][:]; n_t_lon = len(t_lons)
t_lats = t_nc.variables['lat'][:]; n_t_lat = len(t_lats)

#Read in time date and convert to a usable structure 
a_temp = t_nc.variables['temp'][:]
t_rtimes = t_nc.variables['time'][:]
t_times = num2date(t_rtimes, units = t_nc.variables['time'].units , calendar = t_nc.variables['time'].calendar)

#Close File 
t_nc.close()

#Indices of overlap between between precip and MSL data
pr_ii = np.where([t in p_times for t in times])[0]

#Create an array contianing the climatological Mean_MSLP Field
mu = a_p.mean(axis = 0)

#Compute array of anomalies from the mean 
a_pc = a_p - mu

#Create a climatological mean of precip data 
p_mean = a_pr.mean(axis = 0)


#Create a climatological mean of temp data 
t_mean = a_temp.mean(axis = 0)


#Task 1.1 - Plot the Climatological Mean_MSLP_Field over northern Althantic and Europe 
Contour_Plot(mu, "Mean MSLP Field (hPa; Jul/Aug, 1900-2010, daily)")

#Task 1.2 - Calculate Data and Covariance Matrix 
x, x_tilda, s, f, m, columns = data_covariance_matrix()

#Task 1.3 - Compute the Eigenvalues and Eigenvectors of the covariance martix 
lam, e = np.linalg.eig(s)

#Task 1.3 - Resize E (Eigenvectors/Loadings) to plot the first 4 PC loadings
alt_e1 = np.resize(e[:,0], [19,49])
alt_e2 = np.resize(e[:,1], [19,49])
alt_e3 = np.resize(e[:,2], [19,49])
alt_e4 = np.resize(e[:,3], [19,49])

#Task 1.3 - Plot countour plots for 4 PC Loadings 
PCA_Contour_plot(alt_e1, "PCA Loading 1 Mean MSLP Field \n (hPa; Jul/Aug, 1900-2010, daily)")
PCA_Contour_plot(alt_e2, "PCA Loading 2 Mean MSLP Field \n (hPa; Jul/Aug, 1900-2010, daily)")
PCA_Contour_plot(alt_e3, "PCA Loading 3 Mean MSLP Field \n (hPa; Jul/Aug, 1900-2010, daily)")
PCA_Contour_plot(alt_e4, "PCA Loading 4 Mean MSLP Field \n (hPa; Jul/Aug, 1900-2010, daily)")

#Task 1.3 - Rescaled countour plot for 4 PCA Loadings - Scaled to Pressure Anomlilies 
f_inverse = np.linalg.inv(f)
u = np.dot(x_tilda , e)
new_e = np.dot(e, f_inverse)
re_scaled = new_e * np.std(u[:,0])
new_alt_e1 = np.resize(re_scaled[:,0], [19,49])
new_alt_e2 = np.resize(re_scaled[:,1], [19,49])
new_alt_e3 = np.resize(re_scaled[:,2], [19,49])
new_alt_e4 = np.resize(re_scaled[:,3], [19,49])

PCA_Contour_plot_scaled(new_alt_e1, " Rescaled PC Loading 1 Mean MSLP Field \n (hPa; Jul/Aug, 1900-2010, daily)")
PCA_Contour_plot_scaled(new_alt_e2, "Rescaled PC Loading 2 Mean MSLP Field \n (hPa; Jul/Aug, 1900-2010, daily)")
PCA_Contour_plot_scaled(new_alt_e3, "Rescaled PC Loading 3 Mean MSLP Field \n (hPa; Jul/Aug, 1900-2010, daily)")
PCA_Contour_plot_scaled(new_alt_e4, "Rescaled PC Loading 4 Mean MSLP Field \n (hPa; Jul/Aug, 1900-2010, daily)")


# Task 1.4.1 - Calculate U containning PCA Scores within the Rows for U  

u = np.dot(x_tilda , e)

#Task 1.4.1 - Plot Time Series of Principal component Scores for the leading PC Mode.
plot_time_series(times, u[:,0], "g", "Time Series of Principial Component Scores \n for the Leading PC (SNAO) Mode", "Time (Years)", "PCA Score")     

#Task 1.4.2 - Plot a Contour plot for the pressure field for a day where all PC scores equal 0
mean = np.resize(m, [19,49])
Contour_Plot(mean, "Pressure Field for a day where all PC Scores are Equal Zero")

#Calculate Standard Deviation of the SNAO Timeseries 
sd = np.std(u[:,0])

#Task 1.4.2 - Plot a Contour plot for the pressure field for a day where SNAO Score equals SD and all other PC scores equal 0
#Scale Eigen Vectors to allow conversion back to Pressure Values
f_inverse = np.linalg.inv(f)
scaled_eigen_vectors = np.dot(e, f_inverse )

#Create co-efficents using first row of U and replace SNAO Mode with SD
coeff_1= np.full([1,931], sd)
scalling_1 = (coeff_1*scaled_eigen_vectors[:,0])
day_1 = np.resize(m+scalling_1,[19,49])
Contour_Plot(day_1, "Pressure Field for a day where SNAO Score is Equal \n to Standard Deviation of SNAO Time Series")

#Task 1.4.2 - Plot a Contour plot for the pressure field for a day where SNAO Score equals -SD and all other PC scores equal 0
coeff_2= np.full([1,931], -sd)
scalling_2 = (coeff_2*scaled_eigen_vectors[:,0])
day_2 = np.resize(m+scalling_2,[19,49])
Contour_Plot(day_2, "Pressure Field for a day where all SNAO Score is Equal \n to minus Standard Deviation of SNAO Time Series")

#Task 1.5.1 - Plot a Scree plot for the First 12 PCA Modes 
screeval = lam[0:12]
s_mean = 0
s_dev = screeval**0.5
sigma_squared = s_dev**2
sigma_squared_total = sum(sigma_squared)


var_eigen = np.zeros(shape = (12,1))
for i in range (len(var_eigen)):
    var_eigen[i] = (screeval[i]/sigma_squared_total)*100
    
cum_sum = np.cumsum(var_eigen)

fig, ax1 = plt.subplots()
ax1.set_xlabel('Principle Component')
ax2 = ax1.twinx()
ax1.set_ylabel('Standard Deviation', color = 'g')
ax1.plot(s_dev, color = 'g')
ax2.set_ylabel('Cumulative Explained Variance (%)', color = 'b')
ax2.plot(cum_sum, color = 'b')
plt.title("Scree Plot Showing Individual and Cumulative Variance")
plt.grid(True)
fig.tight_layout()
plt.show()

#Task 1.5.2
#Sampling uncertainty under asymtopic assumptions
#Gaussian
n = 6682
mean_lam = np.mean(lam)
sample_est_upper = np.zeros(len(screeval))
sample_est_lower = np.zeros(len(screeval))
#95th Percentile value = 1.96
for j in range(len(screeval)):
    sample_est_lower[j]=(((-1.96/(n/2)**0.5)*screeval[j])+screeval[j])
    
for k in range(len(screeval)):
    sample_est_upper[k]=(((1.96/(n/2)**0.5)*screeval[k])+screeval[k])

sampling_princ = np.zeros((12,3))
sampling_princ[:,0]=np.arange(0,12)
sampling_princ[:,1]=sample_est_lower
sampling_princ[:,2]=sample_est_upper

plt.figure()
plt.plot(screeval, color = "g")
plt.xlabel("Principle Component")
plt.ylabel("Eigenvalue")
plt.title("Scree Plot showing Eigen values for first 12 PC modes \n and theor 95% confidence interval values")
for l in range(len(screeval)):
    plt.vlines(x = sampling_princ[[l],0], ymin = sampling_princ[[l],1], ymax = sampling_princ[[l],2])
    del[l]
plt.vlines(0,0,0, color= "purple", label = "95th Confidence Intervals")
plt.legend()
   
#Task 1.6 - Locate the array index that corresponds to the 20th of July 2007
d_flood = dt.datetime(2007, 7, 20, 12)
fl_i = np.where(times == d_flood)[0][0]

#Plot of Pressure Distribution on the 20th of July 2007 over British Isles
#Showing Depression over a Sepression over south east England and over the Channel 
depression = a_p[6653,:, :]
Plot_Depression(depression, "Pressure Distribution on the 20th of July 2007 over British Isles")

#Set Coefficents using the correct array index from u
#Calculate pressure Field Using Coefficents and scaled Eigen Vectors
#Plot the New Pressure Field over British Isles.
truncated_Depression_plot(5, "Pressure Distribution field approximation \n for the 20th of July 2007 over British Isles \n using 5 Loadings")
truncated_Depression_plot(20, "Pressure Distribution field approximation \n for the 20th of July 2007 over British Isles \n using 20 Loadings")
truncated_Depression_plot(200, "Pressure Distribution field approximation \n for the 20th of July 2007 over British Isles \n using 200 Loadings")
   
# Task 1.7.1 - create composite Plot for temperature and Temperature for pos and neg SNAO

precip_mean_diff, precip_positve_mean_plot, precip_negative_mean_plot, precip_pos, precip_neg = composite_analysis(a_pr)

composition_plot("Precipitation (mm) composite for all days \n when SNAO score is positive", precip_positve_mean_plot, 0, 5, .1)
composition_plot("Precipitation (mm) composite for all days \n when SNAO score is negative", precip_negative_mean_plot, 0, 5, .1)
composition_plot("Precipitation (mm) composite for the difference \n in means between positive and negative SNAO", precip_mean_diff, -1, 1, 0.01)

temp_mean_diff, temp_positve_mean_plot, temp_negative_mean_plot, temp_pos, temp_neg = composite_analysis(a_temp)

composition_plot("Temperature (Degress Celcius) composite for all days \n when SNAO score is positive",temp_positve_mean_plot, 10, 30, 1)
composition_plot("Temperature (Degress Celcius) composite for all days \n when SNAO score is negative", temp_negative_mean_plot, 10, 30, 1)
composition_plot("Temperature (Degress Celcius) composite for the difference \n in means between positive and negative SNAO", temp_mean_diff, -0.7, 1, 0.01)

pr_ttest = np.empty((2736,2))
pr_ttest[0], pr_ttest[1]= stats.stats.ttest_ind(precip_pos[:], precip_neg[:], equal_var = True, nan_policy = 'propagate')
pr_pvalue = pr_ttest[:,1]

# =============================================================================
# pr_stippling = np.zeros((2736))
# for i in range((2736)):
#     if pr_pvalue[i] < 0.05:
#         pr_stippling[i] = pr_pvalue[i]
#     else:
#         pr_stippling[i] = 'nan'
#         
# plot_pr_stippling = np.resize(pr_stippling, [36,76])
# 
# stipple_composition_plot(precip_mean_diff, stippling_variable, plot_pr_stippling, "title")
# stipple_composition_plot(precip_mean_diff, )
# 
# 
# =============================================================================
