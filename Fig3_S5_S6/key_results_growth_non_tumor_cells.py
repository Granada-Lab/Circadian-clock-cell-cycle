#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:53:44 2024

@author: nicagutu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import math

def calculate_stats(series_list):
    combined_values = np.concatenate([series.to_numpy() for series in series_list])  # Flatten the list of Series
    mean = combined_values.mean()
    variance = np.mean((combined_values - mean) ** 2)
    std_dev = math.sqrt(variance)
    return std_dev

plt.rcParams.update({'font.size': 22})    
plt.rcParams['svg.fonttype'] = 'none'

path = 'Data/Population/'
output = 'Figures/Supplementary_figures/Growth_other_cell_lines/'
cellline = 'HEK293T' # NIH3T3, HaCaT, HEK293T
file = f'Circadian_synch_growth_{cellline}.xlsx'
data = pd.read_excel(str(path)+str(file), skiprows=[0])

time = data['Elapsed']
data_cells = data.drop(data.columns[0:2], axis=1)
print(data_cells)

groups = ['A']
labels = {'A': 'high_density'}

# Across inhibitor dose
data = pd.DataFrame(index=data_cells.index.values)
std = pd.DataFrame(index=data_cells.index.values)

for dens in groups:
    df_10 = data_cells[[f'{dens}1', f'{dens}4']]
    df_10 = df_10/df_10.iloc[0]
    data[f'{dens}_inh10'] = savgol_filter(df_10.mean(axis=1), window_length=41, polyorder=5)
    std[f'{dens}_inh10'] = savgol_filter(df_10.std(axis=1), window_length=41, polyorder=5)

    df_5 = data_cells[[f'{dens}2', f'{dens}5']]
    df_5 = df_5/df_5.iloc[0]
    data[f'{dens}_inh5'] = savgol_filter(df_5.mean(axis=1), window_length=41, polyorder=5)
    std[f'{dens}_inh5'] = savgol_filter(df_5.std(axis=1), window_length=41, polyorder=5)

    df_untr = data_cells[[f'{dens}3', f'{dens}6']]
    df_untr = df_untr/df_untr.iloc[0]
    data[f'{dens}_untr'] = savgol_filter(df_untr.mean(axis=1), window_length=41, polyorder=5)
    std[f'{dens}_untr'] = savgol_filter(df_untr.std(axis=1), window_length=41, polyorder=5)
 
for dens in groups:
    fig = plt.figure(figsize=(10,8))
    plt.plot(time, data[f'{dens}_inh10'], label=str('10uM'), linewidth=5)
    lower_bound = data[f'{dens}_inh10'] - std[f'{dens}_inh10']
    upper_bound = data[f'{dens}_inh10'] + std[f'{dens}_inh10']
    plt.fill_between(time, lower_bound, upper_bound, color='tab:blue', alpha=0.3)
    
    plt.plot(time, data[f'{dens}_inh5'], label=str('5uM'), linewidth=5)
    lower_bound = data[f'{dens}_inh5'] - std[f'{dens}_inh5']
    upper_bound = data[f'{dens}_inh5'] + std[f'{dens}_inh5']
    plt.fill_between(time, lower_bound, upper_bound, color='tab:orange', alpha=0.3)
    
    plt.plot(time, data[f'{dens}_untr'], label=str('untreated'), linewidth=5)
    lower_bound = data[f'{dens}_untr'] - std[f'{dens}_untr']
    upper_bound = data[f'{dens}_untr'] + std[f'{dens}_untr']
    plt.fill_between(time, lower_bound, upper_bound, color='tab:green', alpha=0.1)
    
    plt.legend(loc='best')
    plt.xlabel('Time [hours]')
    plt.ylabel('Normalized confluence')
    plt.title(str(cellline)+' '+str(labels[dens]))
    # plt.savefig(output+'Change_dose_'+str(labels[dens])+'_'+str(cellline)+'.svg')
    plt.show()
  


