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
import random
import matplotlib.pyplot as plt # Data Visualization 
import pandas as pd # data processing

# Description Length of a cluster of data
    # DL = 20 + [ nPoints * log2 (sqrt (variance) ) ]
    # distanceArr = the list of distances, where distance = from point to cutting plane
    # start = the first index of distanceArr in the original list
    # end = the last index of distanceArr in the original list
    # 20 = coming from arbitrary resolution (0.01) 

def DL(distanceArr):
    numPoints = len(distanceArr) # number of elements of the distanceArr
    if (numPoints == 0):
        return 0  

    variance = np.var(distanceArr)   
    if (variance == 0): 
        return 20
    elif (math.isnan(variance)):
        return 0
    else:
        return 20 + numPoints * (math.log2(math.sqrt(variance))) ;

def slope(x1, y1, x2, y2):
    return (y2-y1)/(x2-x1)
   
def divide_cluster_DL(dataArr, distanceArr, eigenvector, chopped_clusters, data_plot, dl_plot, show_dl=False, show_hyperplane=False):
    temp_cluster1Dist = distanceArr.copy() 
    temp_cluster2Dist = []
    cut_point = 0
    current_position = 0
    minDL = DL(distanceArr)
       
    dl_arr = [minDL]
    for item in distanceArr:
        temp_cluster1Dist.remove(item) 
        temp_cluster2Dist.append(item)
        current_position = current_position + 1
    
        newDL = (DL(temp_cluster1Dist) 
                 + DL(temp_cluster2Dist))
        
        dl_arr.append(newDL)
        if (newDL < minDL):
            minDL = newDL
            cut_point = current_position
    
    if (show_dl):
        dl_plot.plot(list(range(len(dl_arr))), dl_arr)
        dl_plot.set_xlabel('Index')
        #dl_plot.set_ylabel('Description Length')

    # if no division, return the original data
    if (cut_point == 0): 
        chopped_clusters.append([dataArr.tolist(), minDL])
        return chopped_clusters
    
    # else, divide the data into 2 clusters at the cut_point
    else:
        cluster1 = dataArr[:cut_point]
        cluster2 = dataArr[cut_point:]
        cluster1Dist = distanceArr[:cut_point]
        cluster2Dist = distanceArr[cut_point:]
        
        #print("cut-point", cut_point)
        #print("Cluster 1", cluster1.tolist())
        #print("Cluster 2", cluster2.tolist())
        
        if (show_hyperplane):
            eigvec = [eigenvector[1], -eigenvector[0]]
            #m = slope(0, 0, eigenvector[1], -eigenvector[0])
           
            data_plot.quiver(*(dataArr[cut_point]), *eigvec, color='blue', scale=2)
             
        return (divide_cluster_DL(cluster1, cluster1Dist, eigenvector, [], data_plot, dl_plot, show_dl, show_hyperplane)
                + divide_cluster_DL(cluster2, cluster2Dist, eigenvector ,[], data_plot, dl_plot, show_dl, show_hyperplane))
 
def sorted_dist_list(data, mean, eigenvector):
    # a) Get the distance from the mean (aka, hyperplane) to the origin
    # The hyperplane starts at the mean
    # origin = the starting point of eigenvectors, which is (0,0)
    dist_hyperplane_origin = np.dot(mean.T, eigenvector) # m transpose • vi 

    # b Compute the list of distance from point d to the mean (aka hyperplane)
    dist_hyperplane_point = dist_hyperplane_origin - np.dot(data, eigenvector);
        
    # c) Sort distance_mean_point
    sorted_index = np.argsort(dist_hyperplane_point)[::-1]
    sorted_dist_hyperplane_point = dist_hyperplane_point[sorted_index]
    sorted_data = data[sorted_index]
    return [sorted_dist_hyperplane_point, sorted_data]
        
def chop(data, sorted_eigv, data_plot, dl_plot, show_hyperplane=False):
    # For each eigenvector vi from the sorted list
    input_data = data
    
    for i in range(0, len(sorted_eigv)):
       
        # enter the loop only if not eigenvector 1
        if (isinstance(input_data[0][0], list)):
            clusters = []
            for idx in range(0, len(input_data)):    
                nparray_data = np.array(input_data[idx])
                mean = np.mean(nparray_data, axis = 0)
                sorted_dist_and_data = sorted_dist_list(nparray_data, mean, sorted_eigv[i])
                sorted_dist_hyperplane_point = sorted_dist_and_data[0] 
                sorted_data = sorted_dist_and_data[1]
            
                if (show_hyperplane and len(sorted_eigv) <= 2):
                    sub_clusters = divide_cluster_DL(sorted_data, sorted_dist_hyperplane_point.tolist(), sorted_eigv[i], [], data_plot, dl_plot, True, True)
                else:
                    sub_clusters = divide_cluster_DL(sorted_data, sorted_dist_hyperplane_point.tolist(), sorted_eigv[i], [], data_plot, dl_plot, True, False)

                clusters = sub_clusters + clusters  
                
            input_data = clusters
            print('The number of clusters', sorted_eigv[i], len(clusters))
         
        else:
            mean = np.mean(input_data, axis = 0)
            sorted_dist_and_data = sorted_dist_list(input_data, mean, sorted_eigv[i])
            sorted_dist_hyperplane_point = sorted_dist_and_data[0]
            sorted_data = sorted_dist_and_data[1]
            
            # d) Start sliding the cutting hyper plane
            # Plotting the hyperplane (only accept 2D or less)
            if (show_hyperplane and len(sorted_eigv) <= 2):
                clusters = divide_cluster_DL(sorted_data, sorted_dist_hyperplane_point.tolist(), sorted_eigv[i], [], data_plot, dl_plot, True, True)
            else:
                clusters = divide_cluster_DL(sorted_data, sorted_dist_hyperplane_point.tolist(), sorted_eigv[i], [], data_plot, dl_plot, True, False)

            input_data = [item[0] for item in clusters]
            print('The number of clusters', sorted_eigv[i], len(input_data))
            
    return clusters
    

def merge_cluster_DL(clusters):
    # we will have (K Choose 2) pairs, where K is the number of clusters in C

    if type(clusters) is np.ndarray:
        result = clusters.tolist()
    else:
        result = clusters
        
    for i in range(len(clusters)):
        for x in range(i + 1, len(clusters)):
            a = [*clusters[i], *clusters[x]]
            temp_merge = [*clusters[i][0], *clusters[x][0]]
        
            # Eigenvectors and eigenvalues of the temp_merge
            U, s, VT = LA.svd(temp_merge) 
            eigenvalues, eigenvectors_T = s, VT.T 
    
            # Sort eigenvectors in the descending order of eigenvalues.
            sorted_eig_index = np.argsort(eigenvalues)[::-1] 
            sorted_eigenvectors = eigenvectors_T[:, sorted_eig_index] 
    
            # Distance List
            mean = np.mean(temp_merge, axis = 0)

            sorted_dist_and_data = sorted_dist_list(np.array(temp_merge), mean, sorted_eigenvectors[0])
            sorted_dist_hyperplane_point = sorted_dist_and_data[0] 
                
            dl_add = clusters[i][1] + clusters[x][1]
            dl_union = DL(sorted_dist_hyperplane_point)

            if (dl_union < dl_add):
                
                # remove c1 from clusters
                # print('c1: ', clusters[i])
                result = [item for item in result if item != clusters[i]]
                  
                # remove c2 from clusters
                # print('c2: ', clusters[x])
                result = [item for item in result if item != clusters[x]]
        
                # print('temp_merge: ', temp_merge)
                result.append([temp_merge, dl_union])
               
                return merge_cluster_DL(result)

    # If no merge happens
    return [item[0] for item in result]
 

def main():        
    dataset = pd.read_csv('../data/Mall_Customers.csv', index_col='CustomerID')
    dataset.drop_duplicates(inplace=True)
    data = dataset.iloc[:, [2, 3]].values # ndarray type, each column represents a variable
    data = np.unique(data, axis=0)
    
    fig, (data_plot, dl_plot, clusters_plot) = plt.subplots(1, 3)
    fig.set_size_inches(18.5, 10.5, forward=True)
    
    #data_plot.set_title('Data set with normalized eigenvectors')
    #data_plot.set_xlabel('Annual Income')
    #data_plot.set_ylabel('Spending Score')
    #data_plot.scatter(data[:,0], data[:,1], s = 10, c = 'pink')
    
    # 1. Find the mean of dataset S, axis = 0 means along the column
    data_mean = np.mean(data, axis = 0)

    # 2. Find Co-variance matrix of S
    cov_matrix = np.cov(data, rowvar = False)
    
    # 3. Find the eigenvalues and eigenvectors of Co-variance matrix
    # eigenvalues, eigenvectors = LA.eig(cov_matrix) 
    U, s, VT = LA.svd(data) #SVD version
    eigenvalues, eigenvectors = s, VT.T 
    print(eigenvectors)
    
    # 4. Sort eigenvectors in the descending order of eigenvalues.
    sorted_eig_index = np.argsort(eigenvalues)[::-1] # [::-1] = list[start : stop : step] = reverse a list
    sorted_eigenvectors = eigenvectors[:, sorted_eig_index] # sorted column vectors
    sorted_eigenvectors_T = sorted_eigenvectors.T
    
    sorted_eigenvectors_T = [[0.7709129685111153,0.636940495636273], [0.636940495636273,-0.7709129685111153]]
    dataset = [53.563714462663874,70.43148917119552,68.45589309803958,65.31650471667976,62.659103675706575,62.000271628797925,59.99459537066508,60.79492488238993,61.47482295051258,60.364343912099855,61.4824314898444,62.43728216958551,68.45694850718931,53.563714462663874]
    dl_plot.plot(list(range(len(dataset))), dataset)
    
    # Plotting the eigenvectors  
    #    origin = [0, 0]
    #    data_plot.plot(*data_mean, 'o', label='mean')
    #    for i in range(0, len(sorted_eigenvectors_T)):
    #        color = (random.random(), random.random(), random.random())       
    #        data_plot.quiver(*origin, *(sorted_eigenvectors_T[i]), color=color, 
    #                         label='eigenvector {}'.format(i), scale=1)

	
	
	#dl_plot.plot(list(range(len(datatest))), datatest)
		
    # 5. Chopping process
#    clusters = chop(data, sorted_eigenvectors_T, clusters_plot, dl_plot, True)
#    final_clusters = merge_cluster_DL(clusters)

#    print(len(final_clusters))
    
#    show_cluster = final_clusters
#    for i in range(len(show_cluster)):
#        color = np.array([random.random(), random.random(), random.random()])
#        color=color.reshape(1,-1)
#        cluster = np.array(show_cluster[i])
#        clusters_plot.scatter(cluster[:,0], cluster[:,1], s = 10, c=color)

    data_plot.legend()
    
    
if __name__ == "__main__":
    main()


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
    # Where C = the chopped collections from CHOP()
    For each pair of clusters (Ci, Cj) in C:
        if DL(Ci and Cj) <  DL(Ci) + DL(Cj):
            remove Ci from C
            remove Cj from C
            add (Ci and Cj) to C
            Repeat
    
    if no merge, return C
                
"""

            
    

