#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 09:20:03 2020

@author: niamhmurray
"""

import matplotlib.pyplot as plt 
import numpy as np
import scipy.stats as sps
import statistics as stat
import scipy.special as scsf
from scipy.optimize import fsolve

def read_in_data(filename):
    """
    readindata reads in numerical data from a txt file and produces a 
    numpy array.

    :param filename: the fully qualified pathname of the data file
    :return: the data from the file
    """
    input_data = np.genfromtxt(filename, skip_header=1)
    return (input_data)

def bernoulli_trial(p,x):
    """
    """
    pmf= (p**x)*((1-p)**(1-x))
    return pmf

def mybinomal (x, theta, N):
    """
    """
    coeff = scsf.comb(N,x)
    p = coeff*theta**x *(1-theta)**(N-x)
    return (p) 

def plot_pmf_binomial(om_1, om_2, theta, lambda_t, N, c1, c2, label_1, label_2, title):
    """
    """
    pmf_b = np.zeros(len(om_1))
    
    for j in range(len(om_1)):
        pmf_b[j] = mybinomal(om_1[j], theta, N)
        del j      

    pmf_p = np.zeros(len(om_2))
    
    for i in range(len(om_2)):
        pmf_p[i] = sps.poisson.pmf(om_2[i], lambda_t)
        del i      
    plt.figure()
    plt.plot(om_1, pmf_b, c= c1, label = label_1)
    plt.plot(om_2, pmf_p, c=c2, label = label_2)
    plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.xlabel('# of sucesses')
    plt.ylabel('Probability of Success')

def plot_pmf_geometric_exponential(om_1, om_2, theta, p, c1, c2, label_1, label_2, title):
    """
    """
    pmf_g = np.zeros(len(om_1))
    
    for j in range(len(om_1)):
        pmf_g[j] = sps.geom.pmf(om_1[j], p )
        del j      

    pmf_e = np.zeros(len(om_2))
    
    for i in range(len(om_2)):
        pmf_e[i] = sps.expon.pdf(om_2[i], theta)
        del i      
    plt.figure()
    plt.scatter(om_1, pmf_g, c= c1, label = label_1)
    plt.plot(om_2, pmf_e, c=c2, label = label_2)
    plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.xlabel('Time until first Success')
    plt.ylabel('Probability of such Time')
    
def plot_pmf_normal(om, om_2, mean, sample_mean, sd, sample_sd, c, label, title ):
    """
    """
    pmf = np.zeros(len(om))
    
    for j in range(len(om)):
        pmf[j] = sps.norm(mean, sd).pdf(om[j])
        del j      

        pmf_s = np.zeros(len(om_2))
    
    for i in range(len(om_2)):
        pmf_s[i] = sps.norm(sample_mean, sample_sd).pdf(om_2[i])
        del i 
        
    plt.figure()
    plt.plot(om, pmf, c= c, label = label)
    plt.axvline(mean, c="b", label = "Population Mean")
    plt.scatter(om_2, pmf_s, 10, "r", label = "Sample Elements")
    plt.scatter(sample_mean, sps.norm(sample_mean, sample_sd).pdf(sample_mean), 60, "orange", label = "Sample Mean")
    plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.xlabel('Daily Max Temperature (degrees C)')
#    plt.ylabel('Probability of Occurance')

def plot_pmf_normal1(om, mean, sd, ci, lower, upper, means, c, label, title ):
    """
    """
    pmf = np.zeros(len(om))
    
    for j in range(len(om)):
        pmf[j] = sps.norm(mean, sd).pdf(om[j])
        del j      
    
        plt.figure()
    plt.plot(om, pmf, c= c, label = label)
    plt.axvline(15, c="b", label = "Population Mean")
    #plt.hlines(y, xmin = lower[i], xmax=upper[i], color = "r", linewidth = 0.2)
    plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.xlabel('Daily Max Temperature (degrees C)')
    plt.ylabel('Probability of Occurance')
    
    
def plot_pmf_normal_confidence_intervals(om, om_2, mean, sample_mean, sd, sample_sd, c, label, title, c68, c90, c95, c997, c681, c901, c951, c9971 ):
    """
    """
    pmf = np.zeros(len(om))
    
    for j in range(len(om)):
        pmf[j] = sps.norm(mean, sd).pdf(om[j])
        del j      

        pmf_s = np.zeros(len(om_2))
    
    for i in range(len(om_2)):
        pmf_s[i] = sps.norm(sample_mean, sample_sd).pdf(om_2[i])
        del i 
        
    plt.figure()
    plt.plot(om, pmf, c= c, label = label)
    plt.axvline(mean, c="b", label = "Population Mean")
    plt.axvline(c68, c="r", label = "68% Confidence Interval")
    plt.axvline(c90, c="c", label = "90% Confidence Interval")
    plt.axvline(c95, c="orange", label = "95% Confidence Interval")
    plt.axvline(c997, c="teal", label = "99.7% Confidence Interval")
    plt.axvline(c681, c="r")
    plt.axvline(c901, c="c")
    plt.axvline(c951, c="orange")
    plt.axvline(c9971, c="teal")
    plt.scatter(om_2, pmf_s, 10, "r", label = "Sample Elements")
    plt.scatter(sample_mean, sps.norm(sample_mean, sample_sd).pdf(sample_mean), 60, "orange", label = "Sample Mean")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid(True)
    plt.title(title)
    plt.xlabel('Daily Max Temperature (degrees C)')
    plt.ylabel('Probability of Occurance')    

def confidence_interval(sample_mean, sample_sd, tvalue, m ):
    lower_bound = sample_mean-(tvalue*(sample_sd/(m**(1/2))))
    upper_bound = sample_mean+(tvalue*(sample_sd/(m**(1/2))))
    return ([lower_bound], [upper_bound])

def ci_100(mean, sd, tvalue, m):
    upper_bounds = np.zeros(len(mean))
    lower_bounds = np.zeros(len(mean))
    
    for i in range(len(mean)):
        lower_bounds[i] = mean[i]-(tvalue*(sd[i]/(m**(1/2))))
    
    for j in range(len(mean)):   
        upper_bounds[j] = mean[j]+(tvalue*(sd[i]/(m**(1/2))))
    
    ci = np.zeros([len(mean),2])
    ci[:,0] = lower_bounds
    ci[:,1] = upper_bounds
    
    return ci

def plot_histogram(array, b, color, title, xlabel, ylabel):
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
    plt.figure()
    plt.hist(array,bins=b,color=color, edgecolor='k')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
 
def mme_parameter_estimation(array, x_bar, n):
        summation = np.zeros(n)    
        for i in range(len(array)):
            summation[i] = (array[i]**2)-(x_bar**2)
        base = (1/n)*(sum(summation))
        beta = base/x_bar
        
        alpha = (x_bar**2)/base
        return ([alpha], [beta])

def ml_gamma(x):
    logx = np.log(x)
    x_bar = np.mean(x)
    logx_bar = np.mean(logx)

    def func_alpha(alpha_dum):
        return scsf.digamma(alpha_dum)-np.log(alpha_dum)-logx_bar+np.log(x_bar)

    first_guess = 1
    alpha = fsolve(func_alpha, first_guess)
    beta = x_bar/alpha    

    return alpha,beta

def plot_pmf_gamma(array, alpha1, alpha2, beta1, beta2, c0, c1, c2,label, label1, label2, title ):
    """
    """
        
    plt.figure()
    plt.hist(array, bins=50, density = True, alpha = 0.7, color=c0, label = label )
    plt.plot(np.arange(0,600), sps.gamma.pdf(np.arange(0,600), alpha1, 0, beta1), c= c1, label = label1)
    plt.plot(np.arange(0,600), sps.gamma.pdf(np.arange(0,600), alpha2, 0, beta2), c=c2, label = label2)
    plt.legend()
    plt.grid(True)
    plt.title(title)
    plt.xlabel('Daily Rainfall amounts (mm)')
    plt.ylabel('Probability of Occurance')
        
    
###################Excersise 1###########################
# =============================================================================
# 
# data = read_in_data('wetdays.txt')
# 
# total_wetdays=sum(data[:,2])
# total_days = sum(data[:,1])
# 
# p_wetday = bernoulli_trial(total_wetdays/total_days,1)
# 
# 
# #Binomial Weekly Calculations
# Binom_Week3 = sps.binom.pmf(3,7,p_wetday)
# Binom_Week4 = sps.binom.cdf(4,7,p_wetday)
# Binom_Week2 = 1- sps.binom.cdf(1,7, p_wetday)
# 
# #Poission Weekly Distribution 
# lambda_week= 7*p_wetday
# 
# Possion_Week3 = sps.poisson.pmf(3, lambda_week)
# Possion_Week4 = sps.poisson.cdf(4, lambda_week)
# Possion_Week2 = 1-sps.poisson.cdf(1, lambda_week)
# 
# #Binomial Yearly Calculations
# Binom_Yearly100 = sps.binom.pmf(100,365,p_wetday)
# Binom_Yearly200 = sps.binom.cdf(200,365,p_wetday)
# Binom_Yearly50 = 1- sps.binom.cdf(49,365, p_wetday)
# 
# #Poission Yearly Distribution 
# lambda_year = 365*p_wetday
# 
# Possion_Yearly100 = sps.poisson.pmf(100, lambda_year)
# Possion_Yealry200 = sps.poisson.cdf(200, lambda_year)
# Possion_Yearly50 = 1-sps.poisson.cdf(49, lambda_year)
# 
# #Setting up sample spaces.
# om_week = np.arange(0,8)
# om_year = np.arange(0,366)
# 
# #Plotting PMF of wet days for binomial.
# plot_pmf_binomial(om_week, om_week, p_wetday, lambda_week, 7, "b", "g", "Binomial", "Possion", "Probability Mass function \n of the number of occurances of wet days in a Week ")
# plot_pmf_binomial(om_year, om_year, p_wetday, lambda_year, 365, "b", "g", "Binomial", "Possion", "Probability Mass function \n of the number of occurances of wet days in a Year")
# 
# #Using the Geometric Distribution 
# geometric10 = sps.geom.pmf(10, p_wetday)
# geometric11 = sps.geom.cdf(10, p_wetday)
# geometric31 = sps.geom.cdf(31, p_wetday)
# 
# #Using Exponential Distribution
# theta = 1/lambda_year
# exponential10 = sps.expon.pdf(10,theta)
# exponential11 = sps.expon.cdf(10, theta)
# exponential31 = sps.expon.cdf(31, theta)
# plot_pmf_geometric_exponential(om_year, om_year, theta, p_wetday, "b", "g", "Geometric", "Exponential", "probability Mass/Density Fuction \n of time before the occurance of the first wet day")
# 
# =============================================================================
#####################Excersise 2######################
# =============================================================================
# #Read in Data
# data2 = read_in_data('raindata.txt')
# 
# #Remove missing Values
# rain = data2[:,1]
# mask1 = rain < 600
# valid_rain = rain[mask1]
# valid_time = data2[:,0][mask1]
# 
# #Exploratory Analysis 
# print (sps.describe(valid_rain))
# median = stat.median(valid_rain)
# mode = stat.mode(valid_rain)
# 
# #Plot of rainfall against Time
# plt.figure()
# plt.plot(valid_time, valid_rain, c="b")
# plt.title("Plot of Rainfall amount at Reading University \n against Time")
# plt.xlabel("Time(years)")
# plt.ylabel("Rainfall amount (mm/day)")
# plt.show()
# 
# #Modelling All Rainfall Values
# plot_histogram(valid_rain, 100, "b", "Histogram of rainfall amounts in Reading Observatory \n taking in account all measurments", "Rain Data (mm)", "Frequency")
# 
# plt.figure()
# plt.boxplot(valid_rain)
# plt.yscale('log')
# plt.title("Box plot for rainfall amounts in Reading Observatory \n taking in account all valid measurements")
# plt.show()
# 
# #Modelling non-zero rainfall 
# mask2 = valid_rain > 0
# nonzero_rain = valid_rain[mask2]
# nonzero_rain_time = valid_time[mask2]
# 
# plot_histogram(nonzero_rain, 100, "g", "Histogram of rainfall amounts in Reading Observatory \n taking in account any non-zero measurments", "Non-Zero Rain Data(mm)", "Frequency")
# 
# plt.figure()
# plt.boxplot(nonzero_rain)
# plt.yscale('log')
# plt.title("Box plot for rainfall amounts in Reading Observatory \n taking in account any non-zero measurments")
# plt.show
# 
# #Parameter Estimation MME
# nonzero_rain_mean = np.mean(nonzero_rain)
# mme_parameter = mme_parameter_estimation(nonzero_rain, nonzero_rain_mean, 7201)
# 
# #Parameter estimation MLE
# ml_parameter = ml_gamma(nonzero_rain)
# 
# #
# x=plot_pmf_gamma(nonzero_rain, mme_parameter[0], ml_parameter[0], mme_parameter[1], ml_parameter[1], "r", "g", "b", "Empirical PDF", "MME", "ML", "PDFs of Non-zero Rainfall \n Amounts in Reading Observatory")
# #####################Excersise 3######################
# =============================================================================
#Store Population mean and standard deviation
population_mean = 15
population_sd = 7.0

#Produce and store a normally distributed sample set
np.random.seed(1)
daily_max_temp = np.random.normal(loc= population_mean, scale = population_sd, size = 20)

#Calculate sample mean and sample standard deviation 
sample_mean = np.mean(daily_max_temp)
sample_sd = np.std(daily_max_temp)

#PMF plot of noral distribution 
plot_pmf_normal(np.arange(-15, 50), daily_max_temp, population_mean, sample_mean, population_sd, sample_sd, "g", "pmf", "Daily Max Temp")

# compute T-Test Values
tvalue_68 = sps.t.interval(.68, 19)
tvalue_90 = sps.t.interval(.9, 19)
tvalue_95 = sps.t.interval(.95, 19)
tvalue_997 = sps.t.interval(.997, 19)

#Compute T-Test Statistic Confidence intervals
confidence_68 = confidence_interval(sample_mean, sample_sd, tvalue_68[1], 20 )
confidence_90 = confidence_interval(sample_mean, sample_sd, tvalue_90[1], 20 )
confidence_95 = confidence_interval(sample_mean, sample_sd, tvalue_95[1], 20 )
confidence_997 = confidence_interval(sample_mean, sample_sd, tvalue_997[1], 20 )

# same graph as previous but adding in 4 confidence intervals
plot_pmf_normal_confidence_intervals(np.arange(-15, 50), daily_max_temp, population_mean, sample_mean, population_sd, sample_sd, "g", "pmf", "Daily Max Temp", confidence_68[1], confidence_90[1], confidence_95[1], confidence_997[1],confidence_68[0], confidence_90[0], confidence_95[0], confidence_997[0])

# Calculate Means and sample Standard Deviation for 100 Different sample groups
max_temp_sample_means = np.zeros(100)
max_temp_sample_sd = np.zeros(100)
for i in range (len(max_temp_sample_means)):
    produce_sample = np.random.normal(loc = population_mean, scale = population_sd, size= 20)
    max_temp_sample_means[i]= np.mean(produce_sample)
    max_temp_sample_sd[i]= np.std(produce_sample)

#Plot Histogram of sample means
plot_histogram(max_temp_sample_means, 25, "g", "Histogram Showing Srequency of Sample Means \n From 100 Different Samples of Max Daily Temperature", "Sample Mean Value", "Frequency" )

#Calculate Standard Deviaition of sample Means
sample_means_sd = np.std(max_temp_sample_means)

#Confidence Intervals 
confidence_intervals = ci_100(max_temp_sample_means, max_temp_sample_sd, tvalue_90[1], 20)

#Illustration of pdf and mean, including condidence intervals 
plot_pmf_normal1(np.arange(-15, 50) , population_mean, population_sd, confidence_intervals, confidence_intervals[:,0], confidence_intervals[:,1], max_temp_sample_means, "g", "PDF", "PDF of Max Temperatures \n and Population mean.")
