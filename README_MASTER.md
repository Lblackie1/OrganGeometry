{\rtf1\ansi\ansicpg1252\cocoartf2639
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #Code used for analyses in Blackie et al \'91The Sex of Organ Geometry\'92\
\
##IntensityAlongGutLength\
Code for semi-automatically extracting the xy coordinates and fluorescence intensity of a centreline along the length of a gut from a confocal image in FIJI. Then extracting the mean and standard deviation of fluorescence intensity along gut length in R.\
\
##ProcrustesShapeAnalysis\
Code to process gut 3D centrelines from microCT images into landmarks for geometric morphometric and PCA analysis and procrustes ANOVA.\
\
##ProximityMeasurements\
Code to measure the distance between 3D organ meshes.\
\
##BnlCurvatureCorrelation\
Code to correlate the intensity of fluorescent reporter along gut length (in this case myrGFP reporting bnl) with curvature along gut length. Intensity of fluorescent reporter along gut length is generated from \'91IntensityAlongGutLength\'92 code. Curvature along gut length is generated from \'91CurvatureTorsionMDSAnalysis\'92 code\
\
##CurvatureTorsionMDSAnalysis\
Code to extract features of gut shape such as curvature, torsion and arc length along gut length and plot means and do the mds analysis.}