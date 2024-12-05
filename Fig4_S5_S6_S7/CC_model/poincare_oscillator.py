#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:16:45 2023

@author: nicagutu
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import tqdm
#  import seaborn as sns
#  from scipy.signal import find_peaks

np.random.seed(42)


def coupling(kappa, kind):
#Parameters
    Nosc = 800
    dt = 0.1
    tf = 150

    gg = 1
    ll = 1
    a0 = 1
    acc = 1

#Coupling parameters
    #  kappa = 0.0001 #extracellular circadian coupling
    eps = 0.005 #intracellular circadian - cell cycle

#Vectors initialization
    t = np.arange(0,tf,step=dt)
    X = np.zeros((Nosc,len(t) + 1))
    XX = np.zeros((Nosc,len(t)))
    Y = np.zeros((Nosc,len(t)))
    YY = np.zeros((Nosc,len(t)))

#Initial position
    #  X[:,0] = np.random.uniform(.0,.3,size=(Nosc))
    #  XX[:,0] = np.random.uniform(-1,1,size=(Nosc))
    #  Y[:,0] = np.random.uniform(-.0,.3,size=(Nosc))
    #  YY[:,0] = np.random.uniform(-1,1,size=(Nosc))
    X[:,0] = np.random.uniform(.0,.1,size=(Nosc))
    XX[:,0] = np.random.uniform(-1,1,size=(Nosc))
    Y[:,0] = np.random.uniform(-.0,.1,size=(Nosc))
    YY[:,0] = np.random.uniform(-1,1,size=(Nosc))

#Period distribution
    #  mu1, sigma1 = 23, np.sqrt(6)
    mu1, sigma1 = 22, 3
    period_circ = np.random.normal(mu1, sigma1, Nosc)#+ stats.expon.rvs(0, 1/.1,Nosc)
    mu2, sigma2 = 26, np.sqrt(6)
    #  mu2, sigma2 = 26, 2.5
    period_cell = np.random.normal(mu2, sigma2, Nosc)#+ stats.expon.rvs(0, 1/.1,Nosc)

#  combined_data = np.concatenate((period_circ, period_cell))
#  bin_edges = np.histogram_bin_edges(combined_data, bins=25)

#  plt.figure(figsize=(12,10))
#  num_bins = 25
#  sns.histplot(period_circ, stat='density', kde=True, bins=bin_edges, color='darkblue', label='Circadian')
#  sns.histplot(period_cell, stat='density', kde=True, bins=bin_edges, color='darkgreen', label='Cell cycle')
#  plt.xlabel('Period')
#  plt.legend(loc='best')
#  plt.show()

# Integration over time and population
    SUMS = np.zeros((2,1500))
    for i in tqdm.tqdm(range(len(t) - 1)):
        sum_x = 0
        sum_y = 0

        for m in range(Nosc):
            middle = -gg*(np.sqrt(X[m,i]**2+Y[m,i]**2)-a0)
            X[m,i+1] = X[m,i]+dt*(middle*X[m,i]-2*np.pi*Y[m,i]/period_circ[m]+sum_x*kappa/(2*Nosc))
            Y[m,i+1] = Y[m,i]+dt*(middle*Y[m,i]+2*np.pi*X[m,i]/period_circ[m]+sum_y*kappa/(2*Nosc))
            
            mid = -ll*(np.sqrt(XX[m,i]**2+YY[m,i]**2)-acc)
            XX[m,i+1] = XX[m,i]+dt*(mid*XX[m,i]-2*np.pi*YY[m,i]/period_cell[m]+eps*(X[m,i]+XX[m,i])/2)
            YY[m,i+1] = YY[m,i]+dt*(mid*YY[m,i]+2*np.pi*XX[m,i]/period_cell[m]+eps*(Y[m,i]+YY[m,i])/2)

            for k in range(Nosc):
                sum_x = sum_x+X[k,i] 
                sum_y = sum_y+Y[k,i] 

            SUMS[0,i] = sum_x
            SUMS[1,i] = sum_y


    fig, ax = plt.subplots(1,1, figsize = (8,6))
    for i in range(800):
        X[i,:] = X[i,:]/max(X[i,:])
        #  X[i,:] = np.roll(X[i,:], - np.random.randint(30, 100,1))
        #  inds, _ = find_peaks(X[i,:])
        #  print(np.diff(inds)*dt, period_cell[i], period_circ[i])
        #  X[i,-1] = period_circ[i]
        #  ax.plot(X[i,:], color = 'black', alpha = .08)



#  ax.plot(SUMS[0,:]/np.max(abs(SUMS[0,:])), color = 'red')
#  np.savetxt('half_coupled.csv', X/np.max(X), delimiter = ',')
#  np.savetxt('coupled.csv', X, delimiter = ',')
#  plt.show()
    np.savetxt(f'output/circ_{kind}.csv', X, delimiter = ',')
    np.savetxt(f'output/circ_dist_{kind}.csv', period_circ, delimiter = ',')

