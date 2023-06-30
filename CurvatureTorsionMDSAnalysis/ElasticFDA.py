# -*- coding: utf-8 -*-
"""
Created on Fri May 28 08:37:22 2021

@author: Salem
"""

import numpy as np
import numpy.linalg as la
import numpy.random as npr
#from numpy.polynomial import polynomial as poly

import matplotlib.pyplot as plt


#---------------------------------------------------------------------------------------------------------
from pathlib import WindowsPath
from pathlib import Path
import os
from os import listdir
from os.path import isfile, join

cwd = Path(os.getcwd())

resultsDir = cwd.parents[0]  / "RunResults" 
 
positionsDir = cwd.parents[0]  / "RunResults" / "Smoothed centerlines-cut"

fitResultDir = cwd.parents[0]  / "RunResults" / "FitResults-cut"

fitResultExpDir = cwd / "FitResults-cut"

resultExpDir = cwd/ "RunResults"
#---------------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
import skfda 

from skfda.preprocessing.dim_reduction.projection import FPCA
from skfda.exploratory.stats import var, mean

from skfda.exploratory.stats import fisher_rao_karcher_mean as elastic_mean
from skfda.preprocessing.registration import FisherRaoElasticRegistration as ElasticRegistration
from skfda.preprocessing.registration import landmark_elastic_registration_warping as WarpingRegistration
from skfda.misc.metrics import fisher_rao_distance as ElasticDistance
#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------



dataNames = np.array(['curvatureS']) #, 'angle', 'height',  'torsionS' 
experimentIndices = [5] #1,2,3,4,



runTypes = ['all', 'VC', 'MC', 'VT', 'MT']
colors = ['orange', 'r', 'tab:purple','g', 'tab:blue']


normalizationNames = ["-", "-nlen-"] 

registrationNames = ["-fixed", "-aligned"]

penalty = 0.01


#================================================================================================================================
# distance matrix based on L2 distance for aligned or fixed curves  
#================================================================================================================================
def getDistMatrix(dataIndx=0, experimentIndx=1, nIndx=0, regIndx=0, runIndx=0):
    '''
    
    
    '''
    
    
    dataName = dataNames[dataIndx]
    
    dataMatrix = fdaData(dataIndx=dataIndx, nIndx=nIndx,  experimentIndx=experimentIndx, 
                         regIndx=regIndx, runIndx=runIndx, plotANDsave=True).data_matrix
    
    #return dataMatrix#ElasticDistance(dataMatrix, dataMatrix)
    if runIndx==0:
        numSamples = len(dataMatrix)
      
        
        distMatrix = np.zeros((numSamples,numSamples))
        
        #print(distMatrix.shape)
        
        for ind1 in range(numSamples):
            
            data1 = dataMatrix[ind1]
        
            for ind2 in range(ind1 + 1, numSamples):
          
                data2 = dataMatrix[ind2]
                
                cost = la.norm(data2 - data1)
                
                distMatrix[ind1, ind2] = cost
                distMatrix[ind2, ind1] = cost
                
        
        ndistMatrix = distMatrix/ np.max(distMatrix)     
        
        
        # multidimensional scaling to embed matrix
        mdsPoints = get_MDS_coords(ndistMatrix, experimentIndx=experimentIndx)
        
        plt.matshow(ndistMatrix)
        plt.show()
    
    
        regName = registrationNames[regIndx] 
        normalizationName = normalizationNames[nIndx]   
            
    
    
        
        print("saving:  " + (dataName + regName + normalizationName + "exp-" + 
                                   str(experimentIndx) + "-distMat" ))
    
        np.savetxt(resultExpDir / (dataName + regName + normalizationName + "exp-" + 
                                   str(experimentIndx) + "-distMat.csv"), distMatrix, delimiter=",")
        
        np.savetxt(resultExpDir / (dataName + regName + normalizationName + "exp-" +
                                   str(experimentIndx) + "-mdsCoords.csv"), mdsPoints, delimiter=",")
            
    
    return 
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 



    


#================================================================================================================================ 
# functional data analysis for a data matrix
#================================================================================================================================
def fdaData(dataIndx, nIndx, regIndx, experimentIndx, runIndx, plotANDsave=True):
    '''
    
    
    '''
    
    dataName = dataNames[dataIndx]
    
    
    controlDir = cwd / ("Experiment" + str(experimentIndx)) / "Control-centerlines"
    treatmentDir = cwd / ("Experiment" + str(experimentIndx)) / "Treatment-centerlines"
    
    
    #names0 = orderNames([f[:-4] for f in listdir(positionsDir) if isfile(join(positionsDir, f))])
    
    
    names1 = orderNames([f[:-4] for f in listdir(controlDir) if isfile(join(controlDir, f))])
    
    
    names2 = orderNames([f[:-4] for f in listdir(treatmentDir) if isfile(join(treatmentDir, f))])
    
    names = np.hstack((names1, names2))
    
    # ['all', 'VC', 'MC', 'VT', 'MT']
    runName = runTypes[runIndx]

    
    if runName == 'all':
        names = np.hstack((names1, names2))
    elif runName == 'VC':
        names =  splitNames(names1)[0]
    elif runName == 'MC':
        names =  splitNames(names1)[1]
    elif runName == 'VT':
        names =  splitNames(names2)[0]        
    elif runName == 'MT':
        names =  splitNames(names2)[1]
        
       
    numSamples = len(names)
    
    numGridPoints = 200;
    
    gridPoints = np.linspace(0.0, 1.0, num=numGridPoints)
    
    dataMatrix = np.zeros((numSamples,numGridPoints))
    
    # for ind in range(len(names0)):
        
    #     folder=fitResultDir
        
    #     name = names0[ind]
    
    #     data = np.genfromtxt(folder / (name + "-" + dataName + ".csv"), delimiter=',')
        
    #     arcLength = np.genfromtxt(folder / (name + "-arcLengths.csv"), delimiter=',')
        
    #     nArcLength = (arcLength - min(arcLength))/(max(arcLength) - min(arcLength))
        
    #     #normalization = [arcLength[-1], 1, 1/max(np.abs(data)), arcLength[-1]][dataIndx]
    #     normalization = 1
        
    #     dataMatrix[ind] = reEvaluateArray(nArcLength[:-1], data*normalization, gridPoints)
        
        
        
    for ind in range(len(names)):
        
        folder=fitResultExpDir
        
        name = names[ind]

        
        arcLength = np.genfromtxt(folder / (name + "-arcLengths.csv"), delimiter=',')
        
        nArcLength = (arcLength - min(arcLength))/(max(arcLength) - min(arcLength))
        
            
        data = np.genfromtxt(folder / (name + "-" + dataName + ".csv"), delimiter=',')
        
        
        
        normName = normalizationNames[nIndx]

        if normName == "-":
            
            normalization = 1
            ndata = data*normalization
            
        elif normName == "-nlen-":
            
            normalization = [arcLength[-1], 1, 1/arcLength[-1], arcLength[-1], 1/arcLength[-1]][dataIndx]
            ndata = data*normalization
            
        elif normName == "-nrad-":

            radii = np.genfromtxt(folder / (name + "-radiusS.csv"), delimiter=',')
            #radii = reEvaluateArray(nArcLength, radii, gridPoints)
            
            normalization = [radii, 1, 1/radii, radii, 1/np.mean(radii)][dataIndx]
            ndata = data*normalization
            
        elif normName == "-nn-":
            ndata = (data - min(data))/(max(data) - min(data))
        #
        
        dataMatrix[ind] = reEvaluateArray(nArcLength, ndata, gridPoints)
 
        
    
    
    # for ind in range(len(names2)):
            
    #     folder=fitResultExpDir
        
    #     name = names2[ind]
    
    #     data = np.genfromtxt(folder / (name + "-" + dataName + ".csv"), delimiter=',')
        
    #     arcLength = np.genfromtxt(folder / (name + "-arcLengths.csv"), delimiter=',')
        
    #     nArcLength = (arcLength - min(arcLength))/(max(arcLength) - min(arcLength))
        
        
    #     normName = normalizationNames[nIndx]

    #     if normName == "-":
            
    #         normalization = 1
    #         ndata = data*normalization
            
    #     elif normName == "-nlen-":
            
    #         normalization = [arcLength[-1], 1, 1/arcLength[-1], arcLength[-1], 1/arcLength[-1]][dataIndx]
    #         ndata = data*normalization
            
    #     elif normName == "-nrad-":

    #         radii = np.genfromtxt(folder / (name + "-radiusS.csv"), delimiter=',')
    #         #radii = reEvaluateArray(nArcLength, radii, gridPoints)
            
    #         normalization = [radii, 1, 1/radii, radii][dataIndx]
    #         ndata = data*normalization
            
    #     elif normName == "-nn-":
    #         ndata = (data - min(data))/(max(data) - min(data))
    #     #
        
        
    #     dataMatrix[len(names1) + ind] = reEvaluateArray(nArcLength, ndata, gridPoints)
    

    
    
    fd = skfda.FDataGrid(
        data_matrix=dataMatrix,
        grid_points=gridPoints,
    )


    
    #return fd
    regName = registrationNames[regIndx] 
    
    if regName == "-fixed":
        fd_align = fd
    elif regName == "-aligned":
        elastic_registration = ElasticRegistration(penalty=penalty)

        fd_align = elastic_registration.fit_transform(fd)
    else:
        print("WARNING: registration type not defined")


    
    if plotANDsave:
        
        # fpca_discretized.components_.plot()
        
        fig = fd_align.mean().plot(label="L2 mean", color='black', linestyle='--' )
        elastic_mean(fd_align).plot(fig=fig, label="Elastic mean", color='black')
        
        
        plt.title("Experiment" + str(experimentIndx))
        
        plt.xlabel("arc length")
        plt.ylabel(dataName)
        
        fig.legend()
        
        plt.show()
    
        # fd_align.plot()
        

        
        vfda = var(fd_align)
        mfda = mean(fd_align)
        
        # grid  = np.array(mfda.grid_points[0])
        
        varis = vfda.data_matrix[0]
        means = mfda.data_matrix[0]
        
        
        
        regName = registrationNames[regIndx] 
        normalizationName = normalizationNames[nIndx]   
        
        
        print("saving:  " + dataName + "-" + runName  + regName + normalizationName + "exp-" + str(experimentIndx))
        
        np.savetxt(resultExpDir / (dataName + "-" + runName + regName + normalizationName + "exp-" + 
                                   str(experimentIndx) + "-means.csv"), means, delimiter=",")
        
        np.savetxt(resultExpDir / (dataName + "-" + runName + regName + normalizationName + "exp-" + 
                                   str(experimentIndx) + "-variance.csv"), varis, delimiter=",")
        
        
        
        # plt.plot(grid, means, label="variation")
        # plt.xlabel("arc length")
        # plt.ylabel("mean")
        # plt.axvline(x=mVal, color='r', linestyle='--')
        # plt.show()
        
        
        # plt.plot(grid, varis, label="variation")
        # plt.xlabel("arc length")
        # plt.ylabel("variance")
        # plt.axvline(x=mVal, color='r', linestyle='--')
        # plt.show()
        
    #####   
    
    #fd_align = fd
    return fd_align#, fd_align#, warping
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 



#================================================================================================================================
# functional data analysis for a data matrix
#================================================================================================================================
def fdaData_warping(dataIndx, nIndx,  experimentIndx, runIndx=0):
    '''
    
    
    '''
    
    dataName = dataNames[dataIndx]
    
    
    controlDir = cwd / ("Experiment" + str(experimentIndx)) / "Control-centerlines"
    treatmentDir = cwd / ("Experiment" + str(experimentIndx)) / "Treatment-centerlines"
    
    
    #names0 = orderNames([f[:-4] for f in listdir(positionsDir) if isfile(join(positionsDir, f))])
    
    
    names1 = orderNames([f[:-4] for f in listdir(controlDir) if isfile(join(controlDir, f))])
    
    
    names2 = orderNames([f[:-4] for f in listdir(treatmentDir) if isfile(join(treatmentDir, f))])
    
    

    
    if runIndx == 0:
        names = np.hstack(( names1, names2))
    elif runIndx ==1:
        names =  names1
    
       
    numSamples = len(names)
    
    numGridPoints = 200;
    
    gridPoints = np.linspace(0.0, 1.0, num=numGridPoints)
    
    dataMatrix = np.zeros((numSamples,numGridPoints))
        
        
        
    for ind in range(len(names1)):
        
        folder=fitResultExpDir
        
        name = names1[ind]
    
        data = np.genfromtxt(folder / (name + "-" + dataName + ".csv"), delimiter=',')
        
        arcLength = np.genfromtxt(folder / (name + "-arcLengths.csv"), delimiter=',')
        
        nArcLength = (arcLength - min(arcLength))/(max(arcLength) - min(arcLength))
        
        if nIndx == 0:
            normalization = 1
            ndata = data*normalization
            
        elif nIndx == 1:
            normalization = [arcLength[-1], 1, 1/arcLength[-1], arcLength[-1], 1/arcLength[-1]][dataIndx]
            ndata = data*normalization
        elif nIndx == 2:
            ndata = (data - min(data))/(max(data) - min(data))
        #
        
        dataMatrix[ind] = reEvaluateArray(nArcLength[:-1], ndata, gridPoints)
 
        
    
    
    for ind in range(len(names2)):
            
        folder=fitResultExpDir
        
        name = names2[ind]
    
        data = np.genfromtxt(folder / (name + "-" + dataName + ".csv"), delimiter=',')
        
        arcLength = np.genfromtxt(folder / (name + "-arcLengths.csv"), delimiter=',')
        
        nArcLength = (arcLength - min(arcLength))/(max(arcLength) - min(arcLength))
        
        
        if nIndx == 0:
            normalization = 1
            ndata = data*normalization
            
        elif nIndx == 1:
            normalization = [arcLength[-1], 1, 1/arcLength[-1], arcLength[-1], 1/arcLength[-1]][dataIndx]
            ndata = data*normalization
        elif nIndx == 2:
            ndata = (data - min(data))/(max(data) - min(data))
        #
        
        dataMatrix[len(names1) + ind] = reEvaluateArray(nArcLength[:-1], ndata, gridPoints)
    

    
    
    

    
    fd = skfda.FDataGrid(
        data_matrix=dataMatrix,
        grid_points=gridPoints,
    )
    
 

    elastic_registration = ElasticRegistration(penalty=penalty)

    fd_align = elastic_registration.fit_transform(fd)


    warping = elastic_registration.warping_
    
    
   # plot warping
    fig = warping.plot()

    plt.title(dataName +  "-Experiment" + str(experimentIndx))    
    
    plt.xlabel("normalized arc-length")
    plt.ylabel("warping")
    plt.show()
    
    
    warpingAligned = elastic_registration.fit_transform(warping)
    
    
    fig = warping.mean().plot(label="L2 mean", color=colors[runIndx], linestyle='--' )
    elastic_mean(warpingAligned).plot(fig=fig, label="Elastic mean", color=colors[runIndx])
    
    # plt.axvline(x=mVal, color='black', linestyle='--')
    
    plt.title(dataName +  "-Experiment" + str(experimentIndx))   
    
    plt.xlabel("arc length")
    plt.ylabel("Warping label")
    
    fig.legend()
    
    plt.show()
    
    
    vfda = var(warpingAligned)
    mfda = mean(warpingAligned)
    
    # grid  = np.array(mfda.grid_points[0])
    
    varis = vfda.data_matrix[0]
    means = mfda.data_matrix[0]
    
    
    
    np.savetxt(resultExpDir / (dataName + "-Exp-" + str(experimentIndx) + "-aligned-means-warping.csv"), means, delimiter=",")
    np.savetxt(resultExpDir / (dataName + "-" + "Exp" + str(experimentIndx) + "-aligned-variance-warping.csv"), varis, delimiter=",")
   
   
   
    return warping
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 


#================================================================================================================================
# sum of squares of differences (L2 distance) between functions compared at the same set of arc length values.  
#================================================================================================================================
def analyzeVariance(dataIndx=0,  nIndx=0, regIndx=0, experimentIndx=1, runIndx1=0, runIndx2=0, plotANDsave=True,
                        color='orange', isNormalize=False):
    '''
    
    
    '''
    
    
    
    # get the first and second data matrices
    
    dataMatrix1 = fdaData(dataIndx=dataIndx, nIndx=nIndx, experimentIndx=experimentIndx,
                             regIndx=regIndx, runIndx=runIndx1, plotANDsave=False).data_matrix
    dataMatrix2 = fdaData(dataIndx=dataIndx, nIndx=nIndx, experimentIndx=experimentIndx,
                             regIndx=regIndx, runIndx=runIndx2, plotANDsave=False).data_matrix
    
    numSamples1 = len(dataMatrix1)
    numSamples2 = len(dataMatrix2)
    
    
    numPoints = dataMatrix1.shape[1]
    
    #curve to store distances for each arc-length value
    distances = np.zeros((numPoints,1))
    
    
    for ind1 in range(numSamples1):
        
        data1 = dataMatrix1[ind1]
    
        for ind2 in range(numSamples2):
      
            data2 = dataMatrix2[ind2]
            
            dist = np.abs(data2 - data1)
            
            distances += dist
    
        
    
    
    nDistances = distances/(numSamples1*numSamples2) #np.max(distances)    
    
    
    if isNormalize:
        averageDist = 0.5*(analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
                                           runIndx1=runIndx1, runIndx2=runIndx1,plotANDsave=False) +
                         analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
                                         runIndx1=runIndx2, runIndx2=runIndx2,plotANDsave=False))
            
        
        nDistances = nDistances/averageDist
        
    
    
    if plotANDsave:
        arcS = np.linspace(0, 1, numPoints)
        
        
        plt.plot(arcS, nDistances, color=color)
        
        runType1 = runTypes[runIndx1]
        runType2 = runTypes[runIndx2]
        
        
        plt.title(runType1 + "_vs_" + runType2)
        plt.xlabel("arc length ")
        plt.ylabel("average distance")
        
        plt.show()
        
        dataName = dataNames[dataIndx]
        normName = normalizationNames[nIndx]
        regName = registrationNames[regIndx]
        
        
        np.savetxt(resultExpDir / (dataName + "-" + runType1 + "-" + runType2 + regName + normName + 
                                    "exp-" + str(experimentIndx) + "-ANVA.csv") , nDistances, delimiter=",")
    
            
    
    return nDistances
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 

def analysisOfVariance(dataIndx=0,  nIndx=0, regIndx=0, experimentIndx=1):
    
    numTypes = len(runTypes)
    
    for ind in range(numTypes):
        
        analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
                        runIndx1=ind, runIndx2=ind, color=colors[ind])
    
    
    
    #averageDistance = analyzeVariance(plotANDsave=False)
    
    for ind1 in np.arange(1,numTypes):
        for ind2 in np.arange(ind1 + 1 ,numTypes):
            analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
                            runIndx1=ind1, runIndx2=ind2, color='black', isNormalize=True)
    

    # nDistances1 = 2*analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
    #                 runIndx1=1, runIndx2=1, plotANDsave=False)/(
    #     analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx,experimentIndx=experimentIndx,  
    #                     runIndx1=1, runIndx2=2, plotANDsave=False) + 
    #     analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
    #                     runIndx1=1, runIndx2=3, plotANDsave=False))
    
    
    # nDistances2 = 2*analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
    #                 runIndx1=2, runIndx2=2, plotANDsave=False)/(
    #     analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
    #                     runIndx1=2, runIndx2=1, plotANDsave=False) + 
    #     analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
    #                     runIndx1=2, runIndx2=3, plotANDsave=False))
    
    
    # nDistances3 = 2*analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
    #                 runIndx1=3, runIndx2=3, plotANDsave=False)/(
    #     analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx, experimentIndx=experimentIndx, 
    #                     runIndx1=3, runIndx2=1, plotANDsave=False) + 
    #     analyzeVariance(dataIndx=dataIndx,  nIndx=nIndx, regIndx=regIndx,experimentIndx=experimentIndx,  
    #                     runIndx1=3, runIndx2=2, plotANDsave=False))
    
    # arcS = np.linspace(0, 1, len(nDistances1))

    
    # plt.plot(arcS, nDistances1, color='red')
    # plt.title("within/across")
    # plt.xlabel("arc length ")
    # plt.ylabel("average distance")
    # plt.show()
    
        
    
    # plt.plot(arcS, nDistances2, color='purple')
    # plt.title("within/across")
    # plt.xlabel("arc length ")
    # plt.ylabel("average distance")
    # plt.show()
    
    
    
    # plt.plot(arcS, nDistances3, color='blue')
    # plt.title("within/across")
    # plt.xlabel("arc length ")
    # plt.ylabel("average distance")
    # plt.show()
    
    
    # dataName = dataNames[dataIndx]
    # normName = normalizationNames[nIndx]
    # regName = registrationNames[regIndx]
    
    
    # np.savetxt(resultExpDir / (dataName + "-F" + regName + normName + "relativeDistances.csv"), nDistances1, delimiter=",")
    # np.savetxt(resultExpDir / (dataName + "-V" + regName + normName + "relativeDistances.csv"), nDistances2, delimiter=",")
    # np.savetxt(resultExpDir / (dataName + "-M" + regName + normName + "relativeDistances.csv"), nDistances3, delimiter=",")
    
    return
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================


#================================================================================================================================
# evaluates the array on a given set of points
#================================================================================================================================
def reEvaluateArray(xArray, yArray, gridArray):
    '''
    
    
    '''
    
    newYarray = np.zeros_like(gridArray)
    
    for ind, grid in enumerate(gridArray):
        
       newYarray[ind] = yArray[find_nearest(xArray, grid) ]
    
 
    
    return newYarray
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 


#================================================================================================================================
# the index of the nearest point in the array to the given value
#================================================================================================================================
def find_nearest(array, value):
    array = np.asarray(array)
    
    return (np.abs(array - value)).argmin()
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================
    

from sklearn.manifold import MDS
from scipy.spatial import distance_matrix
#=========================================================================================================================
# Applies multidimensional scaling to visualize the distance matrix
#=========================================================================================================================  
def get_MDS_coords(dist_mat, n_components=3, experimentIndx=1):
    
    
    controlDir = cwd / ("Experiment" + str(experimentIndx)) / "Control-centerlines"
    treatmentDir = cwd / ("Experiment" + str(experimentIndx)) / "Treatment-centerlines"
    
    
    model = MDS(n_components=n_components, dissimilarity='precomputed', metric=True, n_init=4, max_iter=100,
                   verbose=0, eps=0.00001, n_jobs=4, random_state=1)
    out = model.fit_transform(dist_mat)   
    
    return out
    
    xs, ys, zs = np.transpose(out)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    
    names0 = orderNames([f[:-4] for f in listdir(positionsDir) if isfile(join(positionsDir, f))])   
    names1 = orderNames([f[:-4] for f in listdir(controlDir) if isfile(join(controlDir, f))])
    names2 = orderNames([f[:-4] for f in listdir(treatmentDir) if isfile(join(treatmentDir, f))])
    
    
    splittedNames0 = splitNames(names0)
    splittedNames1 = splitNames(names1)
    splittedNames2 = splitNames(names2)
    
    numFemales = len(splittedNames0[0])
    numVirgins = len(splittedNames0[1])
    numMales = len(splittedNames0[2])
    
    numVirCont = len(splittedNames1[1])
    numMalesCont = len(splittedNames1[2])
    
    numVirTreat = len(splittedNames2[1])
    numMalesTreat = len(splittedNames2[2])
    
    colors = np.hstack((np.ones(numFemales), 0.5*np.ones(numVirgins), 0.2*np.ones(numMales), 
                        0.7*np.ones(numVirCont),0.7*np.ones(numMalesCont), 0.0*np.ones(numVirTreat), 
                          0.0*np.ones(numMalesTreat)))
    
    #colors = np.ones(len(xs))
    
    
    ax.scatter(xs, ys, zs, c=colors)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()
    
    return out
#=========================================================================================================================
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#========================================================================================================================= 


#================================================================================================================================
# sum of squares of differences (L2 distance) between functions compared at the same set of arc length values.  
#================================================================================================================================
def orderNames(oldNames): 
    newNames = []

    for name in oldNames:
        if name[0] == "F":
            newNames.append(name)
            
    for name in oldNames:
        if name[0] == "V":
            newNames.append(name)    
        
    for name in oldNames:
        if name[0] == "M":
            newNames.append(name)
        
    return newNames
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 
    

#================================================================================================================================
# sum of squares of differences (L2 distance) between functions compared at the same set of arc length values.  
#================================================================================================================================
def splitNames(oldNames): 
    
    name1 = []
    name2 = []


            
    for name in oldNames:
        if name[0] == "V":
            name1.append(name)    
        
    for name in oldNames:
        if name[0] == "M":
            name2.append(name)
        
    return name1, name2
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 
    



#================================================================================================================================
# sum of squares of differences (L2 distance) between functions compared at the same set of arc length values.  
#================================================================================================================================
def splitLandmarks(oldLandmarks): 
    
    all_names = orderNames([f[:-4] for f in listdir(positionsDir) if isfile(join(positionsDir, f))])
    
    splittedNames = splitNames(all_names)

    numFemales = len(splittedNames[0])
    numVirgins = len(splittedNames[1])
    #numMales = len(splittedNames[2])

 
    landmarks1 = oldLandmarks[:numFemales]
    landmarks2 = oldLandmarks[numFemales:numFemales + numVirgins]
    landmarks3 = oldLandmarks[numFemales + numVirgins:]
        
    
    # print(all_names[:numFemales])
    # print(all_names[numFemales:numFemales + numVirgins])
    # print(all_names[numFemales + numVirgins:])
    
    
    return landmarks1, landmarks2, landmarks3
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 


from scipy.signal import butter, lfilter, freqz
#================================================================================================================================
# making a low-pass filter
#================================================================================================================================
def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 



for i in range(len(dataNames)):
    for j in experimentIndices:
        for k in range(len(normalizationNames)):
            for l in range(len(registrationNames)):
                for n in range(len(runTypes)):
                    getDistMatrix(dataIndx=i, experimentIndx=j, nIndx=k, regIndx=l, runIndx=n)
















