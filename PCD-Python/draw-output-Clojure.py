#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 11:33:51 2021

@author: tphuong
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt # Data Visualization 
import pandas as pd # data processing
import numpy as np
import random

############
### DATA ###

# 1. Full
#original_data = pd.read_csv('../data/Mall_Customers.csv')
full_data = pd.read_csv('../Finding/output_data/original.csv')
full_data = full_data.iloc[:, [0, 1]].values

# 2. Chop
chop_data = pd.read_csv('../Finding/output_data/chop.csv')
# Convert data type string to int
for index, row in chop_data.iterrows():
    if (row['x'] != 'x'):
        row['x'] = float(row['x'])
        row['y'] = float(row['y'])      
chop_data = chop_data.iloc[:, [0, 1]].values

# 3. Merge
merge_data = pd.read_csv('../Finding/output_data/merge.csv')
# Convert data type string to int
for index, row in merge_data.iterrows():
    if (row['x'] != 'x'):
        row['x'] = float(row['x'])
        row['y'] = float(row['y'])      
merge_data = merge_data.iloc[:, [0, 1]].values 
        
###############
### DRAWING ###

# 1. Draw the whole dataset	
fig, (data_plot, chop_plot, merge_plot) = plt.subplots(1, 3)
fig.set_size_inches(20, 10.5, forward=True)
data_plot.scatter(full_data[:,0], full_data[:,1], s = 20, c = 'pink')

# 2. Draw the Merge data	
colors =['brown', 'orange', 'deeppink', 
         'tan','darkseagreen', 'dodgerblue', 
         'rosybrown', 'dimgray', 'teal', 
         'mediumturquoise', 'greenyellow', 'steelblue',
         'thistle', 'blueviolet', 'blue']
start = 0
idxC = 0
merge_arr = merge_data[:, 0]
for i in range(len(merge_arr)):
    if (merge_arr[i] == 'x'):
        merge_plot.scatter(merge_data[start:i, 0], merge_data[start:i,1], s =30, c=colors[idxC])
        #print(merge_data[start:i, 0])
        start = i + 1
        idxC += 1
    
merge_plot.scatter(merge_data[start:, 0], merge_data[start:,1], s = 30, c=colors[idxC])

# 3. Draw the Chop	data
start = 0
chop_arr = chop_data[:, 0]
for i in range(len(chop_arr)):
    if (chop_arr[i] == 'x'):
        color = np.array([random.random(), random.random(), random.random()])
        color = color.reshape(1,-1)
        chop_plot.scatter(chop_data[start:i, 0], chop_data[start:i,1], s = 30, c=color)
        start = i + 1
    
chop_plot.scatter(chop_data[start:, 0], chop_data[start:,1], s = 30, c=color)

	


	
	