#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 12:00:10 2024

@author: nicagutu
"""

import numpy as np
import pandas as pd
from scipy.spatial import Voronoi, voronoi_plot_2d, distance
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.signal import savgol_filter
import seaborn as sns

plt.rcParams.update({'font.size': 22})

# Function to compute areas of Voronoi cells
def compute_voronoi_areas(vor):
    def polygon_area(vertices):
        if len(vertices) < 3:  # Not a polygon
            return 0
        hull = ConvexHull(vertices)
        return hull.volume

    areas = []
    for region in vor.regions:
        if not -1 in region and len(region) > 0:  # Exclude infinite regions
            polygon = [vor.vertices[i] for i in region]
            areas.append(polygon_area(polygon))
    
    return np.array(areas)

# Mean, Variance, Standard Deviation, Coefficient of Variation (CV)
def compute_voronoi_statistics(voronoi_areas):
    mean_area = np.mean(voronoi_areas)
    median_area = np.median(voronoi_areas)
    variance_area = np.var(voronoi_areas)
    sd_area = np.std(voronoi_areas)
    cv_area = (sd_area / mean_area) * 100  # Coefficient of Variation in percentage
    
    # Create a dictionary with all the statistics
    stats = {
        "Statistic": ["Mean", "Median", "Variance", "Standard Deviation", "Coefficient of Variation (CV)"],
        "Value": [mean_area, median_area, variance_area, sd_area, cv_area]
    }
    
    stats_table = pd.DataFrame(stats)
    return stats_table
    
# Function to remove outliers using the IQR method
def remove_outliers(data):
    q1 = np.percentile(data, 25)  # First quartile (25th percentile)
    q3 = np.percentile(data, 75)  # Third quartile (75th percentile)
    iqr = q3 - q1                 # Interquartile range (IQR)
    lower_bound = q1 - 1.5 * iqr   # Lower bound
    upper_bound = q3 + 1.5 * iqr   # Upper bound
    
    # Filter the data within the bounds
    filtered_data = [x for x in data if lower_bound <= x <= upper_bound]
    return filtered_data

path = 'RawData/'
output = 'Figures/'
list_cond = ['100', '50', '25']

statistics_table_all = {}
areas_all = {}
for cond in list_cond:
    centroid_col_df = pd.read_excel(path+'Absolute_centroid_col_'+str(cond)+'.xlsx', index_col=0)
    centroid_row_df = pd.read_excel(path+'Absolute_centroid_row_'+str(cond)+'.xlsx', index_col=0)
       
    # Extract the x and y coordinates, removing the first index column
    # x_coords = centroid_col_df.iloc[0, :].values#.flatten()
    # y_coords = centroid_row_df.iloc[0, :].values#.flatten()
    
    x_coords = centroid_col_df.median(axis=0)
    y_coords = centroid_row_df.median(axis=0)
    
    # plt.plot(y_coords, x_coords, 'o')
    # plt.show()
    
    # Combine x and y coordinates into an array of points
    points = np.column_stack((y_coords, x_coords))
    
    # Compute the Voronoi diagram
    vor = Voronoi(points)
    
    # Plot the Voronoi diagram
    fig, ax = plt.subplots(figsize=(12, 12))
    voronoi_plot_2d(vor, ax=ax, show_vertices=False, line_colors='tab:purple', alpha=0.5, line_width=1.5, point_size=2)
    
    # Display the plot
    plt.title('Voronoi Diagram for '+str(cond)+'% density')
    plt.savefig(output+'Vornoi_diagrams_'+str(cond)+'.svg')
    plt.show()
    
    
    # Compute Voronoi cell areas
    voronoi_areas = compute_voronoi_areas(vor)
    
    # Compute statistics of Voronoi cell areas
    statistics_table = compute_voronoi_statistics(voronoi_areas)
    statistics_table_all[cond] = statistics_table
       
    # Convert areas to numpy array for easier manipulation
    areas = np.array(voronoi_areas)
    
    # Store the areas in the dictionary, keyed by condition
    areas_all[cond] = list(remove_outliers(voronoi_areas))  

    # Plot histogram of areas for this condition
    plt.figure(figsize=(8, 6))
    plt.hist(remove_outliers(voronoi_areas), bins=30, color='skyblue', edgecolor='black', density=True)
    plt.title(f'Distribution of Voronoi Cell Area for {cond}% density')
    # plt.xlim([0,30])
    plt.xlabel('Area of Voronoi Cells')
    plt.ylabel('Density')
    plt.savefig(output+'Distribution_Voronoi_cell_area_'+str(cond)+'.svg')
    plt.show()

print(statistics_table_all)

# Find the maximum length of all the lists
max_length = max(len(areas) for areas in areas_all.values())

# Pad each list with None values to make them the same length
for cond in areas_all:
    areas_all[cond] = areas_all[cond] + [None] * (max_length - len(areas_all[cond]))

# Convert areas_all dictionary into a DataFrame
df = pd.DataFrame(areas_all)

# Plot boxplot and violin plot using seaborn
fig, axs = plt.subplots(1, 2, figsize=(12, 6))

# Boxplot
sns.boxplot(data=df, ax=axs[0])
axs[0].set_title('Boxplot')

# Violin plot
sns.violinplot(data=df, ax=axs[1])
axs[1].set_title('Violin Plot')

# Show the plots
plt.tight_layout()
plt.show()

# Calculate the mean and standard deviation for each condition
means = df.mean()
std_devs = df.std()

# Create a bar plot with error bars showing standard deviation
plt.figure(figsize=(10, 8))
plt.bar(means.index, means.values, yerr=std_devs.values, capsize=5, color='skyblue', edgecolor='black')
plt.xlabel('Seeding density')
plt.ylabel('Mean Voronoi Cell Area [um]')
plt.savefig(output+'Barplot_mean_std_seeding_density_Voronoi_cell_area.svg')
plt.show()

# Convert areas_all dictionary into a long-format DataFrame
df_long = pd.DataFrame([
    {'Condition': cond, 'Area': area}
    for cond, areas in areas_all.items()
    for area in areas if area is not None
])

# Create a strip plot with error bars showing the mean and standard deviation
plt.figure(figsize=(10, 8))
sns.stripplot(x='Condition', y='Area', data=df_long, jitter=True, color='skyblue', alpha=0.5)
sns.pointplot(x='Condition', y='Area', data=df_long, ci='sd', capsize=0.2, color='tab:blue', join=False)
plt.xlabel('Seeding density [%]')
plt.ylabel('Voronoi Cell Area [um^2]')
plt.savefig(output+'Striplot_mean_std_seeding_density_Voronoi_cell_area.svg')
plt.show()


