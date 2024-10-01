# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 11:52:35 2021

@author: gutunn
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size':24})
plt.rcParams['svg.fonttype'] = 'none'

def exponential_func(t, N0, r):
    K = ydata_smoothed[-1]
    return K/(1+((K-N0)/N0)*np.exp(-r*t))

def linear_func(t, a, b):
    return a*t+b


table_path = 'Input/'
output = 'Figures/Supplementary_figures/'

# condition_combine_dict = {101:'high_density_untreated',81:'high_density_0uM',61:'high_density_1.25uM',
#              41:'high_density_2.5uM',21:'high_density_5uM',1:'high_density_10uM',
#               121:'medium_density_untreated',141:'medium_density_0uM',161:'medium_density_1.25uM',
#               181:'medium_density_2.5uM',201:'medium_density_5uM',221:'medium_density_10uM',
#               341:'low_density_untreated',321:'low_density_0uM',301:'low_density_1.25uM',
#               281:'low_density_2.5uM',261:'low_density_5uM',241:'low_density_10uM'}

# condition_combine_dict = {101:'high_density_untreated', 21:'high_density_5uM',1:'high_density_10uM'}
condition_combine_dict = {101:'high_density_untreated', 141:'medium_density_untreated',1:'high_density_10uM'}

# count cell numbers for each condition over frames, and save into a dictionary
count_dict = {}
for j in condition_combine_dict:    
    count_list = [0 for x in range(232)]
    for i in range(20):
        df = pd.read_csv('%s/xy%03d/xy%03d-t_tracking_table.csv' % (table_path, j+i, j+i))
        
        count_list_pos = []
        for k in range(232):
            c = len(set(df.loc[df['frame']==k]['trackId']))
            count_list_pos.append(c)
            
        count_list = [(l+m) for l,m in zip(count_list,count_list_pos)]
    
    count_list = [x/count_list[0] for x in count_list] # this line defines how you normalize cell counts
    count_dict[j] = count_list


# plot results
density = 1
dose = 3

pos_list = [ x for x in condition_combine_dict]
color_list = ['gold','tab:green','tab:blue']
label_list = ['untreated','5uM','10uM'] #'0uM','1.25uM','2.5uM',
label_list = ['high density','50% density','25% density']

print(pos_list)

fit_params1 = []
cov_error1 = []
doubling_times1 = []

fig, ax = plt.subplots()
for j in range(len(label_list)):
    index = pos_list[j]
    ydata = count_dict[index]

    time = np.arange(0,116,0.5)   
    ydata_smoothed = savgol_filter(ydata, window_length=5, polyorder=3)

    tf = np.argmax(ydata_smoothed[0:150])
    time_fit = time[10:tf]
    ydata_fit = ydata_smoothed[10:tf]
    # params, covariance = curve_fit(exponential_func, time_fit, ydata_fit)
    params, covariance = curve_fit(linear_func, time_fit, ydata_fit)

    errors = np.sqrt(np.diag(covariance))
    # y_fit = exponential_func(time, *params)
    y_fit = linear_func(time_fit, *params)
    
    ax.plot(time, ydata_smoothed, label = label_list[j], color = color_list[j],linewidth=5)
    ax.plot(time_fit, y_fit, '--', label='fit: {:.4f} $\pm$ {:.4f}'.format(params[0], errors[0]), alpha=0.5, color = color_list[j],linewidth=5)

    var = 0 # 0 for linear fit and 1 for exp fit
    fit_params1.append(params[var])
    cov_error1.append(errors[var])
    doubling_times1.append(np.log(2)/params[var])

ax.axvline(x=74, color='grey', linestyle='--', alpha=0.5)
ax.legend()
ax.set_xlabel('Time [hours]')
ax.set_ylabel('Normalized number of cells')
fig.set_figheight(12)
fig.set_figwidth(14)
plt.savefig(output+'Cell_counts_changing_density.svg')
plt.show()

fit_params1 = np.array(fit_params1)
cov_error1 = np.array(cov_error1)

p0 = fit_params1[0]
normalized_params = fit_params1 / p0
normalized_errors = normalized_params * np.sqrt((cov_error1 / fit_params1) ** 2 + (cov_error1[0] / p0) ** 2)

plt.figure(figsize=(10,8))
plt.errorbar(label_list, normalized_params, yerr=normalized_errors, fmt='o')
plt.xlabel('Inhibitor concentration') #Inhibitor concentration
plt.ylabel('Growth rate')
# plt.ylim([0.5, 1.1])
# plt.savefig(output+'Logistic_fit_growth_curves_inhibitor.svg')
plt.savefig(output+'Linear_fit_until_saturation_obj_count_density.svg')
plt.show()

print(doubling_times1) 
    