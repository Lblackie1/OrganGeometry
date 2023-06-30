# Code used for analyses in Blackie et al 'The Sex of Organ Geometry'

## IntensityAlongGutLength
Code for semi-automatically extracting the xy coordinates and fluorescence intensity of a centreline along the length of a gut from a confocal image in FIJI. Then extracting the mean and standard deviation of fluorescence intensity along gut length in R.

## ProcrustesShapeAnalysis
Code to process gut 3D centrelines from microCT images into landmarks for geometric morphometric and PCA analysis and procrustes ANOVA.

## ProximityMeasurements
Code to measure the distance between 3D organ meshes.

## BnlCurvatureCorrelation
Code to correlate the intensity of fluorescent reporter along gut length (in this case myrGFP reporting bnl) with curvature along gut length. Intensity of fluorescent reporter along gut length is generated from 'IntensityAlongGutLength' code. Curvature along gut length is generated from 'CurvatureTorsionMDSAnalysis'

## CurvatureTorsionMDSAnalysis
Code to extract features of gut shape such as curvature, torsion and arc length along gut length and plot means and do the mds analysis.
