#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 17:22:01 2020

@author: niamhmurray
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sps
from scipy.special import comb

def read_in_data(filename):
    """
    readindata reads in numerical data from a txt file and produces a 
    numpy array.

    :param filename: the fully qualified pathname of the data file
    :return: the data from the file
    """
    input_data = np.genfromtxt(filename, skip_header=1)
    return (input_data)

def mean_std(array):
    """ produces"""
    x= np.mean(array)
    y= np.std(array)
    return [(x),(y)]

def produce_anomaly(array, mean):
    """produces"""
    anomaly = (array-mean)
    return (anomaly)

def plot_histogram(a, array, b, color, title, xlabel, ylabel):
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
    plt.figure(a)
    plt.hist(array,bins=b,color=color, edgecolor='k')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axvline(array.mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(array.mean()*1.1, max_ylim*0.9, \
             'Mean: {:.2f}'.format(array.mean()))
    
    
def plot_histogram_limited_range(a, array, b, r, c, d, color, title, xlabel, ylabel):
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
    plt.figure(a)
    plt.hist(array,bins=b, range=(c,d), color=color, edgecolor='k')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axvline(array.mean(), color='k', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(array.mean()*1.1, max_ylim*0.9, \
             'Mean: {:.2f}'.format(array.mean()))

def probability_of_value (array, x, N):
    """Produces"""
    count =len(np.where(array == x ) [0])
    prob = count/N
    return (prob)

def probability_percipitation (a, N):
    """calculates"""
    count =len(np.where(a != 0 ) [0])
    prob = (count)/N
    return (prob)

def probability_freezing (b, N):
    """calculates"""
    count =len(np.where(b <0 ) [0])
    prob = (count)/N
    return (prob)

def probaility_of_intersection(N):
    """computes"""
    count =len(np.where((white_christmas[:,2]!=0) & (white_christmas[:,3]<0))[0])
    prob = (count)/N
    return (prob)

def mybinomal (x, theta, N):
    """
    """
    coeff = comb(N,x)
    p = coeff*theta**x *(1-theta)**(N-x)
    return (p) 

def plot_pmf_cmf(om, theta, N, c):
    """
    """
    pmf = np.zeros(len(om))
    cmf = np.zeros(len(om))
    
    for j in range(len(om)):
        pmf[j] = mybinomal(om[j], theta, N)
        cmf[j] = cmf[j-1]+pmf[j]
        del j 
     
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(om, pmf, c= c)
    plt.ylim(-0.1, 1.1)
    plt.xlabel('# of sucesses')
    plt.ylabel('PMF')
    plt.title('binomal dist, N='+str(N)+' \n and theta=' +str(theta))
    
    plt.subplot(1,2,2)
    plt.scatter(om, cmf, c= c)
    plt.ylim(-0.1, 1.1)
    plt.xlabel('# of sucesses')
    plt.ylabel('CMF')
    plt.title('binomal dist, N='+str(N)+' \n and theta=' +str(theta))
    return pmf
    
    
def plot_pmf_cmf_agreements(om, N, pmf_x):
    """
    """
    if (N % 2) == 0:
        Nnew = int(N/2)
        count = [Nnew] 
        
        pmfY = np.zeros(Nnew + 1)
        cmfY = np.zeros(Nnew + 1)
        
        pmfY[0] = pmf_x[Nnew]
        cmfY[0] = pmfY[0]
        
        for j in range (1, Nnew + 1):
            pmfY[j] = 2 * pmf_x[Nnew + j]
            cmfY[j] = cmfY[j - 1] + pmfY[j]
            count.append(Nnew + j)
    
    else: 
        Nnew = int((N + 1)/2)
        count = [Nnew]
    
        pmfY = np.zeros(Nnew)
        cmfY = np.zeros(Nnew)
    
        pmfY[0] = 2 * pmf_x[Nnew]
        cmfY[0] = pmfY[0]
    
        for j in range(1,Nnew):
            pmfY[j] = 2 * pmf_x[Nnew + j]
            cmfY[j] = cmfY[j - 1] + pmfY[j]
            count.append(Nnew + j)
               
         
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(om, pmfY, c= "b")
    plt.ylim(-0.1, 1.1)
    plt.xlabel('# of sucesses')
    plt.ylabel('PMF')
    plt.title('PMF of Agreements \n N='+str(N))
    
    plt.subplot(1,2,2)
    plt.scatter(om, cmfY, c= "b")
    plt.ylim(-0.1,1.1)
    plt.xlabel('# of sucesses')
    plt.ylabel('CMF')
    plt.title("CMF of Agreements \n N="+str(N))
    plt.show()
            
def calculate_pmf_cmf_agreements(N, pmf_x):
    """
    """
    if (N % 2) == 0:
        Nnew = int(N/2)
        count = [Nnew] 
        
        pmfY = np.zeros(Nnew + 1)
        cmfY = np.zeros(Nnew + 1)
        
        pmfY[0] = pmf_x[Nnew]
        cmfY[0] = pmfY[0]
        
        for j in range (1, Nnew + 1):
            pmfY[j] = 2 * pmf_x[Nnew + j]
            cmfY[j] = cmfY[j - 1] + pmfY[j]
            count.append(Nnew + j)
    
    else: 
        Nnew = int((N + 1)/2)
        count = [Nnew]
    
        pmfY = np.zeros(Nnew)
        cmfY = np.zeros(Nnew)
    
        pmfY[0] = 2 * pmf_x[Nnew]
        cmfY[0] = pmfY[0]
    
        for j in range(1,Nnew):
            pmfY[j] = 2 * pmf_x[Nnew + j]
            cmfY[j] = cmfY[j - 1] + pmfY[j]
            count.append(Nnew + j)  
        
        return pmfY, cmfY, count
    
def survival_function(cmfY, om, N):
    """
    """
    survival = np.zeros(len(cmfY))
    survival[0] = cmfY[len(cmfY)-1]
    for j in range (len(cmfY)-1): 
        survival[j + 1] = 1 - cmfY[j]
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(om, survival)
    plt.ylim(-0.1, 1.1)
    plt.xlabel('# of sucesses')
    plt.ylabel('Survival Probablility')
    plt.title('Survival Probability of Agreements \n N='+str(N))
    
    return survival, 
            

X3 = plot_pmf_cmf([0,1,2,3], 0.5, 3, "r")
X4 = plot_pmf_cmf([0,1,2,3,4], 0.5, 4, "r")
X7 = plot_pmf_cmf([0,1,2,3,4,5,6,7], 0.5, 7, "g")
X29 = plot_pmf_cmf([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28, 29], 0.5, 29, "g")

plot_pmf_cmf_agreements([2,3], 3, X3)
plot_pmf_cmf_agreements([2,3,4], 4, X4) 

Y7 = calculate_pmf_cmf_agreements(7, X7)
Y29 = calculate_pmf_cmf_agreements(29, X29)

plot_pmf_cmf_agreements(Y7[2],7, X7)
plot_pmf_cmf_agreements(Y29[2],29, X29)

S7 = survival_function(Y7[1], Y7[2], 7)
S29 = survival_function(Y29[1], Y29[2], 29)

white_christmas= read_in_data('christmas.txt')
rainfal_mean = mean_std(white_christmas[:,2])
mintemp_mean = mean_std(white_christmas[:,3])

rainfall_anomaly = produce_anomaly(white_christmas[:,2], 2.918)
mintemp_anomaly = produce_anomaly(white_christmas[:,3], 2.316)


plot_histogram(10, white_christmas[:,3], 20, "g", "Histogram of Minimum Temp Frequency", "Minimum Temp (Degrees Celcius)", "Frequency")
plot_histogram(11, white_christmas[:,2], 50, "b", "Histogram of Rainfall recorded on Christmas Day", "Rainfall (mm)", "Frequency")

# =============================================================================
# plot_histogram_limited_range(5, white_christmas[:,2], 50, 0, 25, "b", "Histogram of  Rainfall", "Rainfall", "Frequency")
# =============================================================================

plt.figure(12)
plt.hist(white_christmas[:,2],bins=50, range=(0,25), color='r', edgecolor='k')
plt.title('Histomgram of rainfall \n with anomalies removed')
plt.xlabel('Rainfall (mm)')
plt.ylabel('Frequency')
plt.show()


rainfall_skew = sps.skew(white_christmas[:,2])
mintemp_skew = sps.skew(white_christmas[:,3])

rainfall_kurtosis = sps.kurtosis(white_christmas[:,2])
mintemp_kurtosis = sps.kurtosis(white_christmas[:,3])

plt.figure(13)
plt.boxplot(white_christmas[:,3])
plt.title('Boxplot of Christmas Day Temperatures')


probability_2mm = probability_of_value(white_christmas[:,2], 2, 473)

prob_percipitation = probability_percipitation(white_christmas[:,2], 473)

prob_freezing_christmas = probability_freezing(white_christmas[:,3], 473)

prob_intersection = probaility_of_intersection(473)

prob_freeze_given_percip = prob_intersection/prob_percipitation
prob_percip_given_freeze = prob_intersection/prob_freezing_christmas
