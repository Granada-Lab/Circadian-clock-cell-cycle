#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 12:40:43 2024

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyboat import WAnalyzer
from scipy.signal import find_peaks

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

def rescale(signal):
    minimum = np.min(signal)
    maximum = np.max(signal)
    return (signal-minimum)/(maximum-minimum)

def linear_slope(vector):
    x_values = np.arange(len(vector))
    y_values = np.array(vector)
    slope, intercept = np.polyfit(x_values, y_values, 1)
    return slope

path = 'Raw_data/'
output = 'Output/Figures/Supplementary_figures/Circadian_Raw_Signals/'
dose = ['untreated']#, '5uM', '10uM'] #'untreated',, '1.25uM','2.5uM',
density = ['high']#,'medium','low']
channel1 = 'circadian'
channel2 = 'cell_cycle'

dt = 0.5 
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

for i in density:
    for j in dose:

        condition = str(i)+'_density_'+str(j)+'_'
        print(condition)
        
        #Circadian data
        data1 = pd.read_csv(path+condition+channel1+r'.csv', index_col=0)
        #Cell-cycle data
        data2 = pd.read_csv(path+condition+channel2+r'.csv', index_col=0)        
        
        count = 0
        
        for column in data1:
            signal = data1[column].dropna()
            norm_signal = wAn.normalize_amplitude(signal, window_size=50)
            
            signal2 = data2[column].dropna()
            norm_signal2 = wAn.normalize_amplitude(signal2, window_size=50)
            
            peaks1,_ = find_peaks(norm_signal, distance=40)
            peaks2,_ = find_peaks(norm_signal2, distance=40)

            wAn.compute_spectrum(norm_signal, do_plot=False)
            rd = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            power = np.mean(rd['power'])

            wAn.compute_spectrum(norm_signal2, do_plot=False)
            rd2 = wAn.get_maxRidge(power_thresh=0,smoothing_wsize=5)
            power2 = np.mean(rd['power'])
            
            if count<10:
                count += 1
                print(column)
                
                # plt.figure(figsize=(12,10))
                # plt.plot(signal.index.values*0.5, signal, linewidth=5, label='Circadian')
                # # plt.plot(norm_signal2.index.values*0.5, rescale(norm_signal2), linewidth=5, label='Cell cycle')
                # plt.xlabel('Time [hours]')
                # plt.ylabel('Circadian Clock Signal [a.u.]')
                # # plt.legend(loc='best')
                # plt.savefig(output+'Circadian_signal_'+str(count)+'_'+str(condition)+'.svg')
                # plt.show()
                
                plt.figure(figsize=(12,10))
                # plt.plot(norm_signal.index.values*0.5, rescale(norm_signal), linewidth=3, alpha=0.5, label='Circadian')
                # plt.plot(norm_signal2.index.values*0.5, rescale(norm_signal2), linewidth=3, alpha=0.5, label='Cell cycle')
                plt.plot(norm_signal.index.values*0.5, signal+signal2, linewidth=5, label='Overlapped')
                plt.xlabel('Time [hours]')
                plt.ylabel('Overlapped Signal [a.u.]')
                # plt.legend(loc='best')
                # plt.savefig(output+'Overlapped_abs_signal_'+str(count)+'.svg')
                plt.show()
                






