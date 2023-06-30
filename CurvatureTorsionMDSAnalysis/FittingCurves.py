# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 07:47:45 2021

@author: Salem

This scipt takes in a curve as a series of (ordered) points and extracts the curvature and torsion at each point
by fitting the local neighborhood of each point up to cubic order in distance.


Methods: 
    getArcLength(vertices): 
        takes in ordered vertices and computes the arc-length distance from the first vertex 

"""

# find neighborhood of each vertex to be used in approximation
# calculate the arc length for each point
# for each point do the fit up to cubic order. 

import numpy as np
import numpy.linalg as la
import numpy.random as npr
from numpy.polynomial import polynomial as poly


import matplotlib.pyplot as plt


#---------------------------------------------------------------------------------------------------------
#from pathlib import WindowsPath
from pathlib import Path
import os
from os import listdir
from os.path import isfile, join

cwd = Path(os.getcwd())
codeDir = cwd.parents[0]
 

# choose 0 for the full gut centerline and 1 for the midgut or cut centerline. 
cutIndx = 1; cut = ["", "-cut"][cutIndx]
 
positionsDir = cwd  / ("Smoothed centerlines" + cut)
fitResultDir = cwd  /  ("FitResults" + cut)
#fitResultDir = codeDir / "RunResults" / "FitResults-cut"
#---------------------------------------------------------------------------------------------------------


#================================================================================================================================
# get all the files in the positions folder and run fits on them
#================================================================================================================================
def all_fit_runs(resolution=0.05, folder=positionsDir):
    '''
    
    
    '''

    
    fileNames = [f for f in listdir(folder) if isfile(join(folder, f))]
    
    
    for fileName in fileNames:
        
        print(fileName)
        
        verts = np.genfromtxt(folder/fileName, delimiter=',')
        
        print(verts.shape)
        
        find_fit(verts, resolution, fileName) #*getArcLength(verts)[-1]

    
    return 
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 



#================================================================================================================================
# take a list of ordered points and find ranges for the neighborhoods used in the approximation
#================================================================================================================================
def find_fit(vertices, regionSize, fileName):
    '''
    
    
    '''
    
    numVerts = len(vertices)
    
    neibRanges, arc_lengths = find_neighborhood(vertices, regionSize)

    #v_index = 50
    
    curvatures = np.zeros(numVerts)
    torsions = np.zeros(numVerts)
    
    
    for v_index in range(numVerts):
        
        x,y,z = np.transpose(vertices[neibRanges[v_index][0]:neibRanges[v_index][1] + 1] - vertices[v_index])
        
        s_param = (arc_lengths[neibRanges[v_index][0]:neibRanges[v_index][1] + 1] - arc_lengths[v_index])#/(
                  #       arc_lengths[neibRanges[v_index][1] + 1] - arc_lengths[neibRanges[v_index][0]])
   
        # the paramters of the fit using Tikhonov
#        cx = polyRegression(s_param, x, regParam)
#        cy = polyRegression(s_param, y, regParam)
#        cz = polyRegression(s_param, z, regParam)
        
        # the parameters of the fit using numpy polyfit
        cx = poly.polyfit(s_param, x, 3)
        cy = poly.polyfit(s_param, y, 3)
        cz =poly.polyfit(s_param, z, 3)
        
        fitp = np.transpose([cx, cy, cz])
   
        curvatures[v_index] = 2*la.norm(fitp[2])
        
        torsions[v_index] = np.dot(np.cross(fitp[1], 2*fitp[2]), 6*fitp[3])/curvatures[v_index]**2
    
    
    mask = (arc_lengths > regionSize) * ( arc_lengths < (arc_lengths[-1] - regionSize))
    
    plt.plot(arc_lengths[mask], curvatures[mask])
    plt.ylabel('curvature')
    plt.xlabel('arc length')
    plt.show()
    
    #kappa_tau = torsions*curvatures/np.sqrt(np.abs(torsions*curvatures))
    
#    plt.plot(arc_lengths, kappa_tau)
#    plt.ylabel('torsion')
#    plt.xlabel('arc length')
#    plt.show()
    
    
    
    np.savetxt(fitResultDir / (fileName[:-4] + "-curvature.csv"), curvatures, delimiter=",")
    np.savetxt(fitResultDir / (fileName[:-4] + "-torsion.csv"), torsions, delimiter=",")
    np.savetxt(fitResultDir / (fileName[:-4] + "-arcLengths.csv"), arc_lengths, delimiter=",")
    
    return curvatures, torsions
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 



#================================================================================================================================
# take a list of ordered points and find ranges for the neighborhoods used in the approximation
#================================================================================================================================
def find_neighborhood(vertices, regionSize):
    '''
    find_neighborhood(vertices, regionSize):
        takes in vertex positions and a region size and returns ranges of indices for each vertex that
        correspond to the its neighborhood
    
    vertices: (d by numVerts) array representing the positions of the numVerts points in d dimensions.
    regionSize: (arc length) distance on either side of the vertex where the neighboors will be chosen for fitting
    
    '''
    #for how many vertices for each neighborhood?
    
    #check the vertices
    
    numVerts = len(vertices)
    
    s_parameter = getArcLength(vertices)
    
    # the neighborhood range (vertex indices) for each vertex
    neibRanges = np.zeros((numVerts, 2),dtype=int)
    
    for index1 in range(numVerts):
        
        for index2 in range(numVerts):
        
            distA = s_parameter[index1] - s_parameter[np.clip(index1 - index2, 0, numVerts - 1)]
            distB = s_parameter[np.clip(index1 + index2, 0, numVerts - 1)] - s_parameter[index1]
            
            if (distA >= regionSize) or (distB >= regionSize):
                neibRanges[index1] = [np.clip(index1 - index2, 0, numVerts - 1), np.clip(index1 + index2, 0, numVerts - 1) ]
                
                break
            
    return neibRanges, s_parameter
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================
    


#================================================================================================================================
# takes in ordered vertices and computes the arc-length distance from the first vertex 
#================================================================================================================================
def getArcLength(vertices):
    '''
    getArcLength(vertices):
        takes in vertex positions and computes the arc-length parameter for each vertex. 
        this is calculated by adding together all the distances between consecutive vertices. 
    
    vertices: (d by numVerts) array representing the positions of the numVerts points in d dimensions.
    
    
    returns: 
        arcLengths: array of length numVerts composed of the arc length parameter starting from the first vertex
    
    '''
    
    numVerts = len(vertices)
    
    # the neighborhood range (vertex indices) for each vertex
    arcLengths = np.zeros(numVerts)
    
    s_parameter = 0
    
    for index in range(1, numVerts):
        
        s_parameter += la.norm(vertices[index] - vertices[index - 1])
        
        arcLengths[index] = s_parameter
            
    return arcLengths
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================  



#================================================================================================================================
# takes 
#================================================================================================================================
def polyRegression(x_data, y_data, regParam):
    '''
    PolyRegression(data):
        takes in an array of data (xi, yi), i = 1,2,... and returns parmeters of a cubic fit. 
    
    data: (2 by numPoints) array representing the positions of the numPoints data points.
    regParam: The weight of the Tikhonov regulatization
    
    
    returns: 
        list of parameters np.array([a,b,c,d]) corresponding to y(x) = a + b x + c x^2 + d x^3
    
    '''
    
    numPoints = len(x_data)
    
    #Tikhonov matrix
    TikhMatrix = np.diag([0, 0, 2, 1, 1, 1])* regParam
    
    x_data = x_data.reshape((numPoints, 1))
    y_data = y_data.reshape((numPoints, 1))
    
    #input feature matrix
    X_mat = np.hstack((x_data**0, x_data**1, x_data**2, x_data**3, x_data**4, x_data**5))
    
    S1 = np.dot(np.transpose(X_mat), X_mat)
    
    S2 = np.dot(np.transpose(TikhMatrix), TikhMatrix)
    
    params = np.dot(la.pinv(S1 + S2), np.dot(np.transpose(X_mat), y_data))
    
    return params.flatten()

#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================  


#================================================================================================================================
# takes 
#================================================================================================================================
def testRegression(data, regParam):
    '''
    PolyRegression(data):
        takes in an array of data (xi, yi), i = 1,2,... and returns parmeters of a cubic fit. 
    
    data: (2 by numPoints) array representing the positions of the numPoints data points.
    regParam: The weight of the Tikhonov regulatization
    
    
    returns: 
        list of parameters np.array([a,b,c,d]) corresponding to y(x) = a + b x + c x^2 + d x^3
    
    '''
    
    data = np.array([[a/100, (a/100) + (a/100)**3 + 0.2*npr.rand()] for a in range(100)])
    
    numPoints = len(data)
    x_data = data[:, 0].reshape((numPoints, 1))
    y_data = data[:, 1].reshape((numPoints, 1))
    
    plt.plot(x_data, y_data)
    
    params = polyRegression(x_data, y_data, 0.1)
    
    newData = np.array([[a/100, np.dot([1, (a/100), (a/100)**2 ,(a/100)**3], params)] for a in range(100)])
    
    plt.plot(newData[:, 0], newData[:, 1])
    
    return params

#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================  








all_fit_runs()




































