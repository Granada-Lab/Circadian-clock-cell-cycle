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
from scipy.signal import find_peaks
from pyboat import WAnalyzer
from pyboat import ensemble_measures as em

plt.rcParams.update({'font.size': 24})
plt.rcParams['svg.fonttype'] = 'none'

output = 'Figures/Fig2/'
output2 = 'Figures/Supplementary_figures/Circ_period_IMT_distributions/'

dt = 0.1 
lowT = 16
highT = 32
periods = np.linspace(lowT, highT, 200)
wAn = WAnalyzer(periods, dt, time_unit_label='hours')

#Parameters
Nosc = 200
dt = 0.1
tf = 500

gg = 1
ll = 1. #works for [0.5, 1.5]
a0 = 1
acc = 1. # works for [0.5, 1.5]

kappa = 0.7 #extracellular circadian coupling
eps1 = 0.01 #intracellular cell cycle - circadian 
eps2 = 0.0 #intracellular circadian - cell cycle 

# Distribution of parameters only if needed, otherwise all the same

a0_mean = 1
gg_mean = 1

acc_mean = 1 #B0
ll_mean = 1 #gamma

a0 = np.random.normal(a0_mean, 0.1, Nosc)
gg = np.random.normal(gg_mean, 0.1, Nosc)
acc = np.random.normal(acc_mean, 0.1, Nosc) # works for [0.5, 1.5]
ll = np.random.normal(ll_mean, 0.1, Nosc) #works for [0.5, 1.5]

# plt.hist(a0, bins=20)
# plt.show()

#Vectors initialization
t = np.arange(0,tf,step=dt)
X = np.zeros((Nosc,len(t)))
XX = np.zeros((Nosc,len(t)))
Y = np.zeros((Nosc,len(t)))
YY = np.zeros((Nosc,len(t)))

#Initial position
X[:,0] = np.random.uniform(-1,1,size=(Nosc))
XX[:,0] = np.random.uniform(-1,1,size=(Nosc)) 
Y[:,0] = np.random.uniform(-1,1,size=(Nosc))
YY[:,0] = np.random.uniform(-1,1,size=(Nosc))

#Period distribution
mu1, sigma1 = 22, np.sqrt(6)
period_circ = np.random.normal(mu1, sigma1, Nosc)
mu2, sigma2 = 26, np.sqrt(6)
period_cell = np.random.normal(mu2, sigma2, Nosc)
        
# Integration over time and population
for i in range(len(t) - 1):
    sum_x = 0
    sum_y = 0
    for k in range(Nosc):
        sum_x = sum_x+X[k,i] 
        sum_y = sum_y+Y[k,i] 
               
    for m in range(Nosc):
        
        middle = -gg[m]*(np.sqrt(X[m,i]**2+Y[m,i]**2)-a0[m])
        X[m,i+1] = X[m,i]+dt*(middle*X[m,i]-2*np.pi*Y[m,i]/period_circ[m]+sum_x*kappa/(2*Nosc)+eps1*(XX[m,i]))
        Y[m,i+1] = Y[m,i]+dt*(middle*Y[m,i]+2*np.pi*X[m,i]/period_circ[m]+sum_y*kappa/(2*Nosc)+eps1*(YY[m,i]))
        
        mid = -ll[m]*(np.sqrt(XX[m,i]**2+YY[m,i]**2)-acc[m])
        XX[m,i+1] = XX[m,i]+dt*(mid*XX[m,i]-2*np.pi*YY[m,i]/period_cell[m]+eps2*(X[m,i]))
        YY[m,i+1] = YY[m,i]+dt*(mid*YY[m,i]+2*np.pi*XX[m,i]/period_cell[m]+eps2*(Y[m,i]))
        

                   
#%% Phase locking patterns ------------------------------------------

output = '/Users/nicagutu/Nextcloud/Manuscripts/Circadian-CellCycle/Figures/Supplementary_figures/diff_param_cellcycle/'

ph_circ = [((np.arctan2(Y[m, n], X[m, n])+np.pi)/(2*np.pi)) for m in range(Nosc) for n in range(len(t))]
ph_cell = [((np.arctan2(YY[m, n], XX[m, n])+np.pi)/(2*np.pi)) for m in range(Nosc) for n in range(len(t))]
hist, xedges, yedges = np.histogram2d(ph_circ, ph_cell, bins=50)

fig = plt.figure(figsize=(12,10))
plt.imshow(hist.T, extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()], origin='lower', cmap='plasma', interpolation='bilinear')
plt.xlabel(r'Circadian phase ($\theta$/2$\pi$)')
plt.ylabel(r'Cell cycle phase ($\theta$/2$\pi$)')
# plt.title('gamma'+str(ll_mean).replace('.', ',')+'_B'+str(acc_mean).replace('.', ','))
plt.title('eps1_'+str(eps1)+' eps2_'+str(eps2))

# plt.savefig(output+'Ph_lock_patterns_ll'+str(ll_mean).replace('.', '')+'_B0'+str(acc_mean).replace('.', '')+'.svg')
# plt.savefig(output+'Ph_lock_patterns_bidirectional.svg')
# plt.savefig(output+'Ph_lock_patterns_small_cc_influence.svg')
plt.savefig(output+'Ph_lock_patterns_small_circ_influence.svg')

plt.show()

# hist2D = plt.hist2d(ph_circ, ph_cell, bins=25, cmap=plt.cm.jet, density=True)
# print(moransI_square(hist2D[0]))




