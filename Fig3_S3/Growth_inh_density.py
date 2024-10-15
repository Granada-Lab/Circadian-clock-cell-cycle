#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 11:07:03 2023

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from pyboat import WAnalyzer
from pyboat import ensemble_measures as em
import os
import seaborn as sns
    
plt.rcParams.update({'font.size': 28})
plt.rcParams['svg.fonttype'] = 'none'    
  
def exponential_func(t, N0, r):
    K = max(ydata_smoothed)
    return K/(1+((K-N0)/N0)*np.exp(-r*t))

def linear_func(t, a, b):
    return a*t+b

def growthcurve(x0, t, k, kl, LC50, Sm, c, SC50, h): #cell count
    return x0*np.exp(t*k*(1-(Sm*c**h)/(SC50**h+c**h))-t*(kl*c**h)/(LC50**h+c**h))

def exp_func(b, x):
    return 0.068*np.exp(b*x)

def exp_func2(b, x):
    return 0.13*np.exp(-b*x)

dt = 1 
lowT = 16
highT = 40
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')
            
#%%Data
path = 'RawData/'
output = 'Figures/'

file = 'Circadian_growth_TGFbeta.xlsx'

data1_confl = pd.read_excel(path+file, sheet_name=0, skiprows=[0], header=[0])
data1_green = pd.read_excel(path+file, sheet_name=1, skiprows=[0], header=[0]) 

data2_confl = pd.read_excel(path+file, sheet_name=4, skiprows=[0], header=[0])
data2_green = pd.read_excel(path+file, sheet_name=5, skiprows=[0], header=[0]) 

#%%Average results for confluence
data1_confl = data1_confl.drop(columns=['Date Time', 'Elapsed'])
data2_confl = data2_confl.drop(columns=['Date Time', 'Elapsed'])
data_confl = (data1_confl+data2_confl)/2
data_confl = data_confl.drop([2])

line_colors = ['tab:green', 'tab:orange', 'tab:blue']
labels1 = ['10uM', '5uM', 'untreated']
labels2 = ['high', 'medium', 'low']

columns1 = ['A1', 'A2', 'A6']
columns2 = ['A6', 'B6', 'C6']

fit_params1 = []
cov_error1 = []

doubling_times1 = []

plt.figure(figsize=(12,10))
for col, label, color in zip(columns1, labels1, line_colors):
    minimum = data_confl[col].min()
    index = data_confl[col].idxmin()
    ydata = data_confl[col][index:]#/minimum
    
    #shorter range of data for linear fit of the first 24h
    time = ydata.index.values#[0:24]
    ydata = ydata#[0:24]

    ydata_smoothed = savgol_filter(ydata, window_length=9, polyorder=2)
    params, covariance = curve_fit(exponential_func, time, ydata_smoothed)
    # params, covariance = curve_fit(linear_func, time, ydata_smoothed)

    errors = np.sqrt(np.diag(covariance))
    y_fit = exponential_func(time, *params)
    # y_fit = linear_func(time, *params)
    
    var = 1
    fit_params1.append(params[var])
    cov_error1.append(errors[var])
    doubling_times1.append(np.log(2)/params[var])
    
    plt.plot(time, ydata_smoothed, label=label, linewidth=5, color=color)
    plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, color=color, alpha=0.5)
    # plt.plot(time, y_fit, '--', label='fit: {:.2f} $\pm$ {:.2f}'.format(params[0], errors[0]), linewidth=5, color=color, alpha=0.5)

plt.ylabel('Confluence (%)')
plt.xlabel('Time [hours]')
# plt.xticks([0,24,48,72,96,120])
plt.legend(loc='best')
# plt.savefig(output+'Confluence_change_density.svg')
plt.show()

fit_params2 = []
cov_error2 = []

doubling_times2 = []

plt.figure(figsize=(12,10))
for col, label, color in zip(columns2, labels2, line_colors):
    minimum = data_confl[col].min()
    index = data_confl[col].idxmin()
    ydata = data_confl[col][index:]#/minimum
    
    #shorter range of data for linear fit of the first 24h
    time = ydata.index.values#[0:24]
    ydata = ydata#[0:24]

    ydata_smoothed = savgol_filter(ydata, window_length=9, polyorder=2)
    params, covariance = curve_fit(exponential_func, time, ydata_smoothed)
    # params, covariance = curve_fit(linear_func, time, ydata_smoothed)

    errors = np.sqrt(np.diag(covariance))
    y_fit = exponential_func(time, *params)
    # y_fit = linear_func(time, *params)

    var = 1
    fit_params2.append(params[var])
    cov_error2.append(errors[var])
    doubling_times2.append(np.log(2)/params[var])
    
    plt.plot(time, ydata_smoothed, label=label, linewidth=5, color=color)
    plt.plot(time, y_fit, '--', label='exp fit: {:.3f}'.format(params[1]), linewidth=5, color=color, alpha=0.5)
    # plt.plot(time, y_fit, '--', label='fit: {:.2f} $\pm$ {:.2f}'.format(params[0], errors[0]), linewidth=5, color=color, alpha=0.5)

plt.ylabel('Confluence (%)')
plt.xlabel('Time [hours]')
# plt.xticks([0,24,48,72,96,120])
plt.legend(loc='best')
# plt.savefig(output+'Confluence_change_density.svg')
plt.show()

fit_params1 = np.array(fit_params1)
cov_error1 = np.array(cov_error1)
print(fit_params1)
print(cov_error1)

p0 = fit_params1[-1]
normalized_params = fit_params1 / p0
normalized_errors = normalized_params * np.sqrt((cov_error1 / fit_params1) ** 2 + (cov_error1[-1] / p0) ** 2)

plt.figure(figsize=(10,8))
plt.errorbar(labels1, normalized_params, yerr=normalized_errors, fmt='o')
plt.xlabel('Inhibitor concentration') #Inhibitor concentration
plt.ylabel('Relative growth rate')
plt.ylim([0.5, 1.1])
plt.savefig(output+'Logistic_fit_growth_curves_inhibitor.svg')
# plt.savefig(output+'Linear_fit_24h_growth_curves_inhibitor.svg')
plt.show()

fit_params2 = np.array(fit_params2)
cov_error2 = np.array(cov_error2)
print(fit_params2)
print(cov_error2)

p0 = fit_params2[0]
normalized_params = fit_params2 / p0
normalized_errors = normalized_params * np.sqrt((cov_error2 / fit_params2) ** 2 + (cov_error2[0] / p0) ** 2)

plt.figure(figsize=(10,8))
plt.errorbar(labels2, normalized_params, yerr=normalized_errors, fmt='o')
plt.xlabel('Seeding density') #Inhibitor concentration
plt.ylabel('Relative growth rate')
plt.ylim([0.5, 1.1])
plt.savefig(output+'Logistic_fit_growth_curves_density.svg')
# plt.savefig(output+'Linear_fit_24h_growth_curves_density.svg')
plt.show()


    
    
    
    



