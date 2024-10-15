#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 13:21:50 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
import seaborn as sns

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

def remove_outliers(df, column):
    Q1 = df[column].quantile(0.25)  # First quartile (25th percentile)
    Q3 = df[column].quantile(0.75)  # Third quartile (75th percentile)
    IQR = Q3 - Q1  # Interquartile range
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return df[column][(df[column] >= lower_bound) & (df[column] <= upper_bound)]

path = 'RawData/'
output = 'Figures/'
dose = ['untreated', '5uM', '10uM'] # '1.25uM','2.5uM',
density = ['high'] #,'medium','low'
channel1 = 'circadian'
channel2 = 'cell_cycle'

dt = 0.5 
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

#Changing the drug concentration
amplitudes_dict = {}
for j in dose:
    condition = str(density[0])+'_density_'+str(j)+'_'
    print(condition)
    
    ## Circadian data
    # data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
    data1 = pd.read_csv(path+condition+channel1+r'.csv', index_col=0)

    amplitude_circadian = []       
    for column in data1:
        signal = data1[column].dropna()
        detrended = wAn.sinc_detrend(signal, T_c = 50)
        # norm_signal = wAn.normalize_amplitude(signal, window_size=50)
        
        wAn.compute_spectrum(detrended, do_plot=False)
        rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
        power = np.mean(rd['power'])
        amplitude_circadian.extend(rd['amplitude'])
    
    amplitudes_dict[j] = amplitude_circadian

df_ampl = pd.DataFrame({key: pd.Series(value) for key, value in amplitudes_dict.items()}) 
print(df_ampl)

cleaned_df_ampl = pd.DataFrame()
for col in df_ampl:
    cleaned_df_ampl[col] = remove_outliers(df_ampl, col)

plt.figure(figsize = (10, 8))
sns.violinplot(data=cleaned_df_ampl, linewidth=2.5)
plt.ylabel('Amplitude [a.u.]')
plt.xlabel('Condition')
plt.savefig(output+'Amplitude_inhibitor.svg')
plt.show()

plt.figure(figsize = (10, 8))
sns.boxplot(data=cleaned_df_ampl, linewidth=2.5)
plt.ylabel('Amplitude [a.u.]')
plt.xlabel('Condition')
plt.show()
        
#%% Changing the seeding density
dose = ['untreated']#, '5uM', '10uM'] # '1.25uM','2.5uM',
density = ['high' ,'medium','low']

amplitudes_dict = {}
for j in density:
    condition = str(j)+'_density_'+str(dose[0])+'_'
    print(condition)
    
    ## Circadian data
    # data1 = pd.read_csv(path+condition+channel1+r'_filtered.csv', index_col=0)
    data1 = pd.read_csv(path+condition+channel1+r'.csv', index_col=0)

    amplitude_circadian = []       
    for column in data1:
        signal = data1[column].dropna()
        detrended = wAn.sinc_detrend(signal, T_c = 50)
        # norm_signal = wAn.normalize_amplitude(signal, window_size=50)
        
        wAn.compute_spectrum(detrended, do_plot=False)
        rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
        power = np.mean(rd['power'])
        amplitude_circadian.extend(rd['amplitude'])
    
    amplitudes_dict[j] = amplitude_circadian

df_ampl = pd.DataFrame({key: pd.Series(value) for key, value in amplitudes_dict.items()}) 

cleaned_df_ampl = pd.DataFrame()
for col in df_ampl:
    cleaned_df_ampl[col] = remove_outliers(df_ampl, col)
    

plt.figure(figsize = (10, 8))
sns.violinplot(data=cleaned_df_ampl, linewidth=2.5)
plt.ylabel('Amplitude [a.u.]')
plt.xlabel('Condition')
plt.savefig(output+'Amplitude_density.svg')
plt.show()

plt.figure(figsize = (10, 8))
sns.boxplot(data=cleaned_df_ampl, linewidth=2.5)
plt.ylabel('Amplitude [a.u.]')
plt.xlabel('Condition')
plt.show()
        
