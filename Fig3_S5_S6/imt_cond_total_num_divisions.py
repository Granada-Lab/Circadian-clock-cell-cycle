#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 12:50:35 2024

@author: nicagutu
"""

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

path = 'RawData/'
output = 'Figures/'

#%% Changing TGF-beta inhibitor dose
dose = ['untreated', '5uM', '10uM']
density = 'high'#['high', 'medium', 'low']

# Create an empty dictionary to store intermitotic times (IMT) for each division group across all doses
imt_all_conditions = {str(num): {dose_level: [] for dose_level in dose} for num in range(2, 6)}

for i in dose:
    condition = f'{density}_density_{i}_'
    print(condition)

    # Load data
    data = pd.read_excel(path + f'divisions_{condition}.xlsx')

    num_divisions = []
    for col in data:
        div_events = data[col]
        division_times = np.where(div_events == 1)[0]
        
        num_divisions.append(len(division_times))
        if 2 <= len(division_times) <= 5:
            imt_single = []
            for ii in range(len(division_times) - 1):
                imt_single.append((division_times[ii + 1] - division_times[ii]) * 0.5)
            # Append the average IMT for this particular cell to the corresponding division group for this dose
            imt_all_conditions[str(len(division_times))][i].append(np.mean(imt_single))

# Now let's plot boxplots per division group, with different conditions in each boxplot
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

for idx, num_div in enumerate(range(2, 6)):
    # Collect the data for the current division group
    division_group_data = [imt_all_conditions[str(num_div)][dose_level] for dose_level in dose]
    
    # Create a boxplot
    axes[idx].boxplot(division_group_data, labels=dose, showfliers=False)
    axes[idx].set_title(f'{num_div} divisions')
    axes[idx].set_ylim([15, 50])
    axes[idx].set_ylabel('IMT [hours]')

plt.tight_layout()
plt.savefig(output+'IMT_dose_condition_total_divisions.svg')
plt.show()    

#%% Changing seeding density

dose = 'untreated'#['untreated', '5uM', '10uM']
density = ['high', 'medium', 'low']

# Create an empty dictionary to store intermitotic times (IMT) for each division group across all doses
imt_all_conditions = {str(num): {dens_level: [] for dens_level in density} for num in range(2, 6)}

for i in density:
    condition = f'{i}_density_{dose}_'
    print(condition)

    # Load data
    data = pd.read_excel(path + f'divisions_{condition}.xlsx')

    num_divisions = []
    for col in data:
        div_events = data[col]
        division_times = np.where(div_events == 1)[0]
        
        num_divisions.append(len(division_times))
        if 2 <= len(division_times) <= 5:
            imt_single = []
            for ii in range(len(division_times) - 1):
                imt_single.append((division_times[ii + 1] - division_times[ii]) * 0.5)
            # Append the average IMT for this particular cell to the corresponding division group for this dose
            imt_all_conditions[str(len(division_times))][i].append(np.mean(imt_single))

# Now let's plot boxplots per division group, with different conditions in each boxplot
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

for idx, num_div in enumerate(range(2, 6)):
    # Collect the data for the current division group
    division_group_data = [imt_all_conditions[str(num_div)][dens_level] for dens_level in density]
    
    # Create a boxplot
    axes[idx].boxplot(division_group_data, labels=density, showfliers=False)
    axes[idx].set_title(f'{num_div} divisions')
    axes[idx].set_ylim([15, 50])
    axes[idx].set_ylabel('IMT [hours]')

plt.tight_layout()
plt.savefig(output+'IMT_density_condition_total_divisions.svg')
plt.show()    
