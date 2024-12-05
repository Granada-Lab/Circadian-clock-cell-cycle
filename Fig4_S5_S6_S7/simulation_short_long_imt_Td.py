#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 18:12:48 2024

@author: nicagutu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:45:04 2023

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pyboat import WAnalyzer

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

output = 'Figures/Supplementary_figures/Circ_period_IMT_distributions/'

dt = 0.1 # the sampling interval, 0.5hours
lowT = 10
highT = 40
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

#Parameters
Nosc = 200
dt = 0.1
tf = 200

gg = 1
ll = 1
a0 = 1
acc = 1

kappa = 0 #extracellular circadian coupling # 0.7 on and 0 off
eps1 = 0. #intracellular circadian - cell cycle #0
eps2 = 0.1 #intracellular circadian - cell cycle #0.01 or 0.015 or 0.2 for stonger effects

kappa_vect = [0, 0.7]
period_mean_cc = [16, 34]

df_periods_cc = pd.DataFrame()
count = 0
for jj in range(len(period_mean_cc)):

    for kk in range(len(kappa_vect)):
        
        kappa = kappa_vect[kk]
        
        #Vectors initialization
        t = np.arange(0,tf,step=dt)
        X = np.zeros((Nosc,len(t)))
        XX = np.zeros((Nosc,len(t)))
        Y = np.zeros((Nosc,len(t)))
        YY = np.zeros((Nosc,len(t)))
        
        #Initial position
        np.random.seed(0)
        X[:,0] = np.random.uniform(-1,1,size=(Nosc))
        XX[:,0] = np.random.uniform(-1,1,size=(Nosc)) #-0.5 if eps=0.1 and -0.3 if eps=0.15
        Y[:,0] = np.random.uniform(-1,1,size=(Nosc))
        YY[:,0] = np.random.uniform(-1,1,size=(Nosc))
        
        #Period distribution
        mu1, sigma1 = 24, np.sqrt(6)
        period_circ = np.random.normal(mu1, sigma1, Nosc)
        mu2, sigma2 = period_mean_cc[jj], np.sqrt(6)
        period_cell = np.random.normal(mu2, sigma2, Nosc)
            
        # Integration over time and population
        for i in range(len(t) - 1):
            sum_x = 0
            sum_y = 0
            
            for k in range(Nosc):
                sum_x = sum_x+X[k,i] 
                sum_y = sum_y+Y[k,i] 
        
            for m in range(Nosc):
                middle = -gg*(np.sqrt(X[m,i]**2+Y[m,i]**2)-a0)
                X[m,i+1] = X[m,i]+dt*(middle*X[m,i]-2*np.pi*Y[m,i]/period_circ[m]+sum_x*kappa/(2*Nosc)+eps1*(XX[m,i]))
                Y[m,i+1] = Y[m,i]+dt*(middle*Y[m,i]+2*np.pi*X[m,i]/period_circ[m]+sum_y*kappa/(2*Nosc)+eps1*(YY[m,i]))
                
                mid = -ll*(np.sqrt(XX[m,i]**2+YY[m,i]**2)-acc)
                XX[m,i+1] = XX[m,i]+dt*(mid*XX[m,i]-2*np.pi*YY[m,i]/period_cell[m]+eps2*(X[m,i]))
                YY[m,i+1] = YY[m,i]+dt*(mid*YY[m,i]+2*np.pi*XX[m,i]/period_cell[m]+eps2*(Y[m,i]))
        
                    
        #%% Circdian and cell cycle oscillations
        phases_circ = []
        periods_circ_end = []
        # plt.figure(figsize=(12,10))
        for ii in range(Nosc):
            signal = np.sin(np.arctan2(Y[ii], X[ii]))
            # plt.plot(t, signal, linewidth=5)
            wAn.compute_spectrum(signal, do_plot=False)
            rd = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)
            periods_circ_end.extend(rd['periods'].values.tolist())  
        # plt.xlabel('Time [hours]')
        # plt.ylabel('Circadian clock oscillations')
        # # plt.savefig(output+'Phases_circ_extr'+str(kappa).replace('.', '')+'_intra'+str(eps).replace('.', '')+'_Nosc'+str(Nosc)+'.svg')
        # plt.show()
        
        periods_cc_end = []
        # plt.figure(figsize=(12,10))
        for ii in range(Nosc):
            signal2 = np.sin(np.arctan2(YY[ii], XX[ii]))
            # plt.plot(t, signal2, linewidth=5)
            wAn.compute_spectrum(signal2, do_plot=False)
            rd2 = wAn.get_maxRidge(power_thresh=0, smoothing_wsize=5)
            periods_cc_end.extend(rd2['periods'].values.tolist())
        
        
        #%% Instantaneous periods
        combined_data = np.concatenate((periods_circ_end, periods_cc_end))
        bin_edges = np.histogram_bin_edges(combined_data, bins=25)
        
        plt.figure(figsize=(12,10))
        num_bins = 25
        sns.histplot(periods_circ_end, stat='density', kde=True, kde_kws={'bw_adjust': 5}, bins=bin_edges, color='darkblue', label='Circadian end')
        sns.histplot(periods_cc_end, stat='density', kde=True, kde_kws={'bw_adjust': 5}, bins=bin_edges, color='darkgreen', label='Cell cycle end')
        plt.xlabel('Period')
        plt.legend(loc='best')
        plt.show() 
        
        print(np.median(periods_cc_end))
        
        df_periods_cc[count] = periods_cc_end
        count += 1

print(df_periods_cc)

means = df_periods_cc.mean()
stds = df_periods_cc.std()

ln2 = np.log(2)
r_mean = ln2 / means
r_error = (ln2 / (means**2)) * stds

print("Mean growth rates:\n", r_mean)
print("Uncertainty in growth rates:\n", r_error)

relative_r_mean = r_mean.copy()
relative_error = r_mean.copy()

reference_indices = [0, 2]  # Reference columns for each group
for ref in reference_indices:
    group_indices = range(ref, ref + 2)  # Columns [0,1] for ref=0 and [2,3] for ref=2
    for i in group_indices:
        relative_r_mean[i] = r_mean[i] / r_mean[ref]
        relative_error[i] = relative_r_mean[i] * np.sqrt(
            (r_error[i] / r_mean[i])**2 + (r_error[ref] / r_mean[ref])**2)

# Display results
print("Relative r_mean:\n", relative_r_mean)
print("Relative errors:\n", relative_error)

df_periods_cc.columns = ['Uncoupled \n IMT<Tcirc', 'Coupled \n IMT<Tcirc', 'Uncoupled \n IMT>Tcirc', 'Coupled \n IMT>Tcirc']

x_positions = np.arange(len(df_periods_cc.columns))

plt.figure(figsize=(12, 10))
plt.errorbar(df_periods_cc.columns, relative_r_mean, yerr=relative_error, fmt='o', color='blue', ecolor='black',
              elinewidth=2, capsize=6)
plt.xlabel('Extracellular circadian coupling')
plt.ylabel('Relative growth rates')
plt.savefig(output+'Simulated_relative_growth_IMT_bigger_smaller_Tcirc.svg')
plt.show()



