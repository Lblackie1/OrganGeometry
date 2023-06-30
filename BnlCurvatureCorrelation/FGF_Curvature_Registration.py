import numpy as np
#import numpy.linalg as la
#import numpy.random as npr
#from numpy.polynomial import polynomial as poly

import matplotlib.pyplot as plt


from scipy.stats import spearmanr, pearsonr


#---------------------------------------------------------------------------------------------------------
#from pathlib import WindowsPath
from pathlib import Path
import os
from os import listdir
from os.path import isfile, join

cwd = Path(os.getcwd())

fgfResultsDir = cwd /  "Results"

positionsDir = cwd  / "Raw Positions"
fgfDataDir = cwd  /  "fgfData"
syntheticDir = cwd  /  "SyntheticData"
fitResultDir = cwd  /  "Curvatures"


dataSources = [fgfDataDir, syntheticDir]
sourceNames = ["fgfData", "synthetic"]
#---------------------------------------------------------------------------------------------------------


import skfda 

from skfda.preprocessing.dim_reduction.projection import FPCA
#from skfda.preprocessing.registration import ElasticRegistration
#from skfda.preprocessing.registration.ElasticRegistration import elastic_mean
from skfda.exploratory.stats import var, mean

from skfda.exploratory.stats import fisher_rao_karcher_mean as elastic_mean
from skfda.preprocessing.registration import FisherRaoElasticRegistration as ElasticRegistration


dataNames = np.array(['curvatureS'])
runTypes = ['V', 'M']


penalty = 0.01

#================================================================================================================================
# sum of squares of differences (L2 distance) between functions compared at the same set of arc length values.  
#================================================================================================================================
def getCorrelations(dataIndx = 0 , runIndx=0,  dataSourceIndx = 0, plotANDsave=False):
    '''
    
    
    '''
    
    dataName = dataNames[dataIndx]
    runType = runTypes[runIndx]
    
    dataDir = dataSources[dataSourceIndx]
    
    sourceName = sourceNames[dataSourceIndx]
    
    
    saveName  = dataName + '-' + sourceName + '-' + runType
    

    #takes all the names ordered (F, V, M)
    if(runType == 'all'):        
        allNames = orderNames([f[:-4] for f in listdir(positionsDir) if isfile(join(positionsDir, f))])
        fgfNames =  orderNames([f[:-4] for f in listdir(dataDir) if isfile(join(dataDir, f))])
        
    #choses one of F, V, M
    elif (runType == 'F'):
        allNames = splitNames([f[:-4] for f in listdir(positionsDir) if isfile(join(positionsDir, f))])[0]
        fgfNames =  splitNames([f[:-4] for f in listdir(dataDir) if isfile(join(dataDir, f))])[0]
    elif (runType == 'V'):
        allNames = splitNames([f[:-4] for f in listdir(positionsDir) if isfile(join(positionsDir, f))])[1]
        fgfNames =  splitNames([f[:-4] for f in listdir(dataDir) if isfile(join(dataDir, f))])[1]
    elif (runType == 'M'):
        allNames = splitNames([f[:-4] for f in listdir(positionsDir) if isfile(join(positionsDir, f))])[2]
        fgfNames =  splitNames([f[:-4] for f in listdir(dataDir) if isfile(join(dataDir, f))])[2]
    
    numSamples = len(allNames)  
    
    numFGFs = len(fgfNames)
    
    #numGridPoints = 200;
    
    #gridPoints = np.linspace(0.0, 1.0, num=numGridPoints)
    
    #dataMatrix = np.zeros((numSamples,numGridPoints))
    
    
    
    #numSamples =10
    sampleRange = np.arange(numSamples)
    #fgfRange = np.arange(13, numFGFs);
    fgfRange = np.arange(numFGFs);
    
    correlationCoeffs = []
    
    for ind1 in sampleRange:
        
        curvatureData = np.genfromtxt(fitResultDir / (allNames[ind1] + "-" + dataName + ".csv"), delimiter=',')
        
        #arcLength = np.genfromtxt(fitResultDir / (name + "-arcLengths.csv"), delimiter=',')
        
        #curvatureData = np.stack((curvatures,arcLength), axis = 1)
       

        for ind2 in fgfRange:
             
     
            fgfData  = np.genfromtxt(dataDir / (fgfNames[ind2] + ".csv"), delimiter=',')
         

            fd = alignTwoData(fgfData, curvatureData, plotANDsave=plotANDsave)
            
            #return fgfData
            
            correlationCoeffs.append(fd)
    
    
    np.savetxt(fgfResultsDir / (saveName + ".csv"), correlationCoeffs, delimiter=",")
        
    #####    
    return np.array(correlationCoeffs)
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 




#================================================================================================================================
#   find the 
#================================================================================================================================
def alignTwoData(data1, data2, numGridPoints=200, plotANDsave=False):
    '''
    
    
    '''
    

    
    gridPoints = np.linspace(0.0, 1.0, num=numGridPoints)
    
    dataMatrix = np.zeros((2,numGridPoints))
    
 
        
        
    x1, y1 = data1.T
        
    nx1 = (x1 - min(x1))/(max(x1) - min(x1)) 
    ny1 = (y1 - min(y1))/(max(y1) - min(y1))  
        
    dataMatrix[0] = reEvaluateArray(nx1, ny1, gridPoints)
    
    #return dataMatrix[0]
    
    
    x2, y2 = data2.T
        
    nx2 = (x2 - min(x2))/(max(x2) - min(x2))    
    ny2 = (y2 - min(y2))/(max(y2) - min(y2))  
        
    dataMatrix[1] = reEvaluateArray(nx2, ny2, gridPoints)
    

    
    fd = skfda.FDataGrid(
        data_matrix=dataMatrix,
        grid_points=gridPoints,
    )
    
    

    g, f = fd[0], fd[1]
    


    elastic_registration = ElasticRegistration(template=g, penalty=penalty)

    
    f_align = elastic_registration.fit_transform(f)
    
    
    corr, _ = pearsonr(f_align.data_matrix.flatten(),  g.data_matrix.flatten())
    
   
    
    if plotANDsave:
        
    
        fig = g.plot(color="orange")
        f_align.plot(fig=fig, color='b')
        #f.plot(fig=fig, color='green', linestyle='dashed')
        
        
        # Legend
        fig.axes[0].legend(['$fgf$', r'$k \circ \gamma $']) #, r'$k$'
        plt.show()
        
    
        
    #####    
    return corr
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



#================================================================================================================================
# sum of squares of differences (L2 distance) between functions compared at the same set of arc length values.  
#================================================================================================================================
def orderNames(oldNames): 
    newNames = []

    for name in oldNames:
        if name[0]== "F":
            newNames.append(name)
            
    for name in oldNames:
        if name[0]== "V":
            newNames.append(name)    
        
    for name in oldNames:
        if name[0]== "M":
            newNames.append(name)
        
    return newNames
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================ 
    

#================================================================================================================================
# sum of squares of differences (L2 distance) between functions compared at the same set of arc length values.  
#================================================================================================================================
def splitNames(oldNames): 
    
    names1 = []
    names2 = []
    names3 = []

    for name in oldNames:
        if name[0]== "F":
            names1.append(name)
            
    for name in oldNames:
        if name[0]== "V":
            names2.append(name)    
        
    for name in oldNames:
        if name[0]== "M":
            names3.append(name)
        
    return names1, names2, names3
#================================================================================================================================
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#================================================================================================================================




for i in range(len(dataSources)):
    for j in range(len(runTypes)):
        for k in range(len(dataNames)):
            getCorrelations(dataSourceIndx=i, runIndx=j, dataIndx=k)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            