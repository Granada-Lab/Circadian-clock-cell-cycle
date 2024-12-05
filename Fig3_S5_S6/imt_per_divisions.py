#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 12:55:32 2024

@author: nicagutu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyboat import WAnalyzer
from pyboat import ensemble_measures as em
from pyboat import plotting as pl 
import seaborn as sns
import random
from scipy.stats import kurtosis
import scipy.stats as stats
from scipy.stats import ttest_ind

plt.rcParams.update({'font.size': 22})    
plt.rcParams['svg.fonttype'] = 'none'

path = '/Users/nicagutu/Nextcloud/Manuscripts/Circadian-CellCycle/Data/1stExp_ilastik/Division_matrix/'
output = '/Users/nicagutu/Nextcloud/Manuscripts/Circadian-CellCycle/Figures/Supplementary_figures/IMT_divisions_TGFb_Density/'

#%% Changing TGF-beta inhibitor dose
dose = ['untreated', '5uM', '10uM']
density = 'high'#['high', 'medium', 'low']

for i in dose:
    condition = f'{density}_density_{i}_'
    print(condition)

    # Load data
    data = pd.read_excel(path + f'divisions_{condition}.xlsx')
    
    num_divisions = []
    imt = {str(num): [] for num in range(2, 6)}
    imt_all = []
    for col in data:
        div_events = data[col]
        division_times = ((np.where(div_events==1)[0]))
        
        num_divisions.append(len(division_times))
        imt_single = []
        if 2 <= len(division_times) <= 5:
            for ii in range(len(division_times)-1):
                imt_single.append((division_times[ii+1]-division_times[ii])*0.5)
            imt[str(len(division_times))].append(np.mean(imt_single))
            # imt_all.append(np.mean(imt_single))
            
    df_imt = pd.DataFrame.from_dict(imt, orient='index')
    df_imt = df_imt.transpose()
    df_imt.columns = [str(num)+' divisions' for num in range(2, 6)]
    
    print('Mean value ', np.mean(df_imt.mean()))
    print('Median value ', np.mean(df_imt.median()))
    print('Coefficient of variation of means ', np.std(df_imt.mean())/np.mean(df_imt.mean()))
    
    # for jj in range(4):
    #     coljj = df_imt.columns[jj]
    #     print('n ',len(df_imt[coljj].dropna()))
    #     for ii in range(4):
    #         colii = df_imt.columns[ii]
    #         if ii > jj:
    #             if jj+1 == ii:
    #                 t_stat, p_value = ttest_ind(df_imt[coljj].dropna(), df_imt[colii].dropna(), equal_var=False)
    #                 print("P-value of "+str(coljj)+" to "+str(colii)+":", p_value)

    # for jj in range(1,4):
    #     coljj = df_imt.columns[jj]
    #     t_stat, p_value = ttest_ind(df_imt[coljj].dropna(), df_imt[df_imt.columns[0]].dropna(), equal_var=False)
    #     print("P-value of "+str(coljj)+" to "+str(df_imt.columns[0])+":", p_value)

    fig = plt.figure(figsize=(10,10))
    df_imt.boxplot(showfliers=False)
    plt.ylim([15,50])
    plt.ylabel('IMT [hours]')
    plt.title(str(condition))
    # plt.savefig(output+'IMT_boxplot_numdivisions_'+str(condition)+'.svg')
    plt.show()
    

#%% Changing seeding density

dose = 'untreated'#['untreated', '5uM', '10uM']
density = ['high', 'medium', 'low']

for i in density:
    condition = f'{i}_density_{dose}_'
    print(condition)

    # Load data
    data = pd.read_excel(path + f'divisions_{condition}.xlsx')
    
    num_divisions = []
    imt = {str(num): [] for num in range(2, 6)}
    imt_all = []
    for col in data:
        div_events = data[col]
        division_times = ((np.where(div_events==1)[0]))
        
        num_divisions.append(len(division_times))
        imt_single = []
        if 2 <= len(division_times) <= 5:
            for ii in range(len(division_times)-1):
                imt_single.append((division_times[ii+1]-division_times[ii])*0.5)
            imt[str(len(division_times))].append(np.mean(imt_single))
            # imt_all.append(np.mean(imt_single))
            
    df_imt = pd.DataFrame.from_dict(imt, orient='index')
    df_imt = df_imt.transpose()
    df_imt.columns = [str(num)+' divisions' for num in range(2, 6)]
    
    print('Mean value ', np.mean(df_imt.mean()))
    print('Median value ', np.mean(df_imt.median()))
    print('Coefficient of variation of means ', np.std(df_imt.mean())/np.mean(df_imt.mean()))
    
    # for jj in range(4):
    #     coljj = df_imt.columns[jj]
    #     print('n ',len(df_imt[coljj].dropna()))
    #     for ii in range(4):
    #         colii = df_imt.columns[ii]
    #         if ii > jj:
    #             if jj+1 == ii:
    #                 t_stat, p_value = ttest_ind(df_imt[coljj].dropna(), df_imt[colii].dropna(), equal_var=False)
    #                 print("P-value of "+str(coljj)+" to "+str(colii)+":", p_value)

    # for jj in range(1,4):
    #     coljj = df_imt.columns[jj]
    #     t_stat, p_value = ttest_ind(df_imt[coljj].dropna(), df_imt[df_imt.columns[0]].dropna(), equal_var=False)
    #     print("P-value of "+str(coljj)+" to "+str(df_imt.columns[0])+":", p_value)

    fig = plt.figure(figsize=(10,10))
    df_imt.boxplot(showfliers=False)
    plt.ylim([15,50])
    plt.ylabel('IMT [hours]')
    plt.title(str(condition))
    # plt.savefig(output+'IMT_boxplot_numdivisions_'+str(condition)+'.svg')
    plt.show()
    

    

    