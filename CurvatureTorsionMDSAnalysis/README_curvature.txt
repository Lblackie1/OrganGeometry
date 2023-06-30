This code requires Mathemtica (tested on Mathemtica 13.1 for windows) and python 3.9 (tested in spyder, anaconda3 for windows)

Python packages used (installed using pip or conda): numpy, scipy, matplotlib, pathlib, os, sklearn, skfda


*** make sure the data is arranged so that the mathemtica and python code are in the same directory as the data folders.

*** For the experiments, each data folder has the subdirectories "Control-centerlines" and "Treatment-centerlines".

*** The data folders have the names Experiment#, where # is the experiment number that may need to be specified in the code.

*** The Mathemtica file "AnalyzingExperiments.nb" should be run first to (a) extract the full gut vertices (results saved in the directory "smoothed centerlines") or the midgut (saved in "smoothed centerlines-cut").

*** Next run the python code "FittingCurves.py" (from the terminal or in spyder you can hit the run file). This will generate a list of quantities such as curvature, torsion, 
arc-length that are saved in the folder "FitResults" (for the full gut analysis) or "FitResults-cut" (for the midgut analysis).

*** Next run the block titled (*Show Curves and results*) in the file "AnalyzingExperiments.nb". This defined several useful functions to visualize and smooth the data. 
For example, tracePlot[][#] will plot the centerline points color coded according to arclength. 

*** Next run the block (*Generate smoothed data*). 

*** Next to perform the mds analysis, run the file "ElasticFDA.py". 

*** Then back in "AnalyzingExperiments.nb" we can plot the mean curvature for control and experiment flies and generate the mds results.


The results for the other experiments and oregonR datasets follow in a similar manner.
 