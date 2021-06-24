#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 11:40:11 2021

@author: tphuong
"""
# PSEUDOCODE for CHOP(S)
import numpy as np # linear algebra
from numpy import linalg as LA
import math
import matplotlib.pyplot as plt # Data Visualization 

import pandas as pd # data processing
dataset = pd.read_csv('./Mall_Customers.csv', index_col='CustomerID')
dataset.drop_duplicates(inplace=True)
X = dataset.iloc[:, [2, 3]].values # ndarray type (subset of data for testing purpose)

plt.scatter(X[:,0], X[:,1], s = 50, c = 'pink', label = 'Original cluster')

# Description Length of a cluster of data
    # DL = 20 + [ nPoints * log2 (sqrt (variance) ) ]
    # distanceArr = the list of distances, where distance = from point to cutting plane
    # start = the first index of distanceArr in the original list
    # end = the last index of distanceArr in the original list
    # 20 = coming from arbitrary resolution (0.01) 

def DL(distanceArr, start, end):
    numPoints = end - start # number of elements of the distanceArr
    variance = np.var(distanceArr)    
    if (variance == 0): 
        return 20
    else:
        return 20 + numPoints * (math.log2(math.sqrt(variance))) ;

        
def divide_cluster_DL(dataArr, distanceArr, chopped_clusters):
    temp_cluster1Dist = distanceArr.copy() 
    temp_cluster2Dist = []
    cut_point = 0
    current_position = 0
    minDL = DL(distanceArr, 0, len(distanceArr))
    
    for item in distanceArr:
        temp_cluster1Dist.remove(item) 
        temp_cluster2Dist.append(item)
        current_position = current_position + 1
        newDL = DL(temp_cluster1Dist, current_position, len(distanceArr)) 
        + DL(temp_cluster2Dist, 0, current_position) 
        
        if (newDL < minDL):
            minDL = newDL
            cut_point = current_position
    
    # if no division, return the original data
    if (cut_point == 0): 
        chopped_clusters.append(dataArr.tolist())
        return chopped_clusters
    
    # else, divide the data into 2 clusters at the cut_point
    else:
        cluster1 = dataArr[:cut_point]
        cluster2 = dataArr[cut_point:]
        cluster1Dist = distanceArr[:cut_point]
        cluster2Dist = distanceArr[cut_point:]
    
        #print('CLUSTER 1: ', cluster1, '\n')
        #print('CLUSTER 2: ', cluster2, '\n')
        return divide_cluster_DL(cluster1, cluster1Dist, []) + divide_cluster_DL(cluster2, cluster2Dist, [])

    
def chop(data):
    # S = a set of data points (array: each column represents a variable)
    # S_mean = the mean of S
    # S_CovMatrix = co-variance matrix of S
    
    # 1. Find the mean of dataset S, axis = 0 means along the column
    data_mean = np.mean(data, axis = 0)
    
    # 2. Find Co-variance matrix of S
    cov_matrix = np.cov(data, rowvar = False)
    
    # 3. Find the eigenvalues and eigenvectors of C
    eigenvalues, eigenvectors = LA.eig(cov_matrix)
    
    # Plotting the eigenvectors
    eig_vec1 = eigenvectors[:,0]
    eig_vec2 = eigenvectors[:,1]    
    plt.quiver(*data_mean, *eig_vec1, color=['r'], scale=10, label='eigenvector 1', alpha=0.5)
    plt.quiver(*data_mean, *eig_vec2, color=['b'], scale=10, label='eigenvector 2' ,alpha=0.5)
    
    
    # 4. Sort both eigenvalues and eigenvectors in the descending order of eigenvalue.
    sorted_eig_index = np.argsort(eigenvalues)[::-1] # [::-1] = list[start : stop : step] = reverse a list
    sorted_eigenvalue = eigenvalues[sorted_eig_index]
    sorted_eigenvectors = eigenvectors[:, sorted_eig_index] # sorted column vectors
        
    #print('original data: \n', data)
    
    # 5. For each eigenvector vi from the sorted list
    for vec in sorted_eigenvectors.T:
        # a) Establish the cutting hyper plane 
        cutting_hyper_plane = np.dot(data_mean.T, vec) # m transpose • vi 

        # b Compute the list of distance from point d to the cutting plane n
        distance_from_plane = cutting_hyper_plane - np.dot(data, vec);

        
        # c) Sort the points in order of distance from the cutting hyper plane
        sorted_index = np.argsort(distance_from_plane)[::-1]
        sorted_distance_from_plane = distance_from_plane[sorted_index]
        sorted_data = data[sorted_index]
      
        
        # d) Start sliding the cutting hyper plane
        # print('DATA', sorted_distance_from_plane, '\n')
        # sorted_distance_from_plane and sorted_data have the same length and index
        clusters = divide_cluster_DL(sorted_data, sorted_distance_from_plane.tolist(), [])   
        print()
        print('eigenvector', vec, len(clusters))
        #print('data points, sorted based on their distances to the hyper plane: ', sorted_data)
        print('clusters', clusters)


def merge_cluster_DL(clusters):
    # we will have (K Choose 2) pairs, where K is the number of clusters in C
    if type(clusters) is np.ndarray:
        clusters = clusters.tolist()
        result = clusters
    else:
        result = clusters
         
    for i in range(len(clusters)):
        for x in range(i + 1, len(clusters)):
            #print(clusters[i], clusters[x])
            temp_merge = [*clusters[i], *clusters[x]]
            
            if (DL(temp_merge) < (DL(clusters[i]) + DL(clusters[x]))):   
                
                #print('c1: ', clusters[i])
                result = [item for item in result if item != clusters[i]]
                #print('Result after deleting c1: ', result)
                  
                #print('c2: ', clusters[x])
                result = [item for item in result if item != clusters[x]]
                #print('Result after deleting c2: ', result)
        
                #print('temp_merge: ', temp_merge)
                result.append(temp_merge)
                #print('Final output: ', result)
                
                return merge_cluster_DL(result)
            
    # If no merge happens
    return clusters
   
"""
Testing  
# 1. Chop
clustersAter = chop(X) 
 
# 2. plot chopped clusters
print('Number of clusters:', len(clustersAter))
for i in clustersAter:
    i = np.array(i)
    rgb = np.random.rand(3,)
    plt.scatter(i[:,0], i[:,1], s = 50, color=[rgb], label = 'Cluster 1', alpha=0.5)  
   
# 3. Merge   
#a = np.array([[1,1], [2,2], [3,3], [4,4]])
#merge_cluster_DL(a)

"""


# PSEUDOCODE for CHOP(S)
"""
# CHOP results in some groups of data points being divided unnecessarily 
def CHOP(S):
    # S = a set of n-dimensional data points
    # m = the mean of S
    # C = co-variance matrix of S
    
    1. Find mean of S
    2. Find Co-variance matrix of S
    3. Find the eigenvalues and eigenvectors of C
    4. Sort both eigenvalues and eigenvectors in the descending order of eigenvalue.
    5. For each eigenvector vi from the sorted list (4)
    
        a) Establish the cutting hyper plane: The cutting plane is perpendicular to the eigenvector vi. To start, we will pick the one that passes through meean m.
            CuttingPlane n = m transpose • vi 
            
    
        b) Compute the list of distance from point d to the cutting plane n
            For any point d, distance = n - d • vi
            
        c) Sort the points in order of distance from the cutting hyper plane
            
        d) Start sliding the cutting hyper plane
            A = the sorted list of data points S
            B = empty list
            cutPoint = 0
            position = 0 # a tracking index
            minDL = DL(A) -- DL(A) is the description length of the entire data set A 
            
            For each point dj in A do:
                - Add dj to B
                - Remove dj from A
                - Increment the position: position = position + 1
                - newDL = DL(A) + DL(B)
                
                if newDL < minDL:
                    minDL = newDL
                    cutPoint = position
            
            if cutPoint > 0:
                - divide A into 2 collections S1 and S2 at the cutPoint index
            
            Recursively apply CHOP to S1 and S2
            
            Return a list of chopped collections: [S1, S2, ... , Sn]
            
        e) If A remains unchange, the data point cannot be represented with a smaller DL 
            return S 
"""    


# PSEUDOCODE for MERGE(C)
"""
# MERGE(C)
# check if the DL would be reduced if they were to be merged (handle non-convex)
def MERGE(C): 
    # C = the chopped collections from CHOP(S)
    For each pair of clusters (Ci, Cj) in C:
        if DL(Ci and Cj) <  DL(Ci) + DL(Cj):
            C = [original C without Ci, Cj and merge(Ci and Cj)]
            Repeat
    
    if no merge, return C
    
"""

            
    
