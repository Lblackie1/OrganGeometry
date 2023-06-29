#Package installation:
#install.packages("geomorph")
#install.packages("nat")
#install.packages("ggfortify")
#install.packages("ggplot2")
#If needed for the XML package:
#install.packages("XML", type = "binary")

library(nat)
library(geomorph)
library(XML)
library("ggfortify")
library("ggplot2")

#Load gut centrelines and convert to landmarks
setwd("~/DemoData_handdcrtraRNAi")
# read gut traces into R as neurons. 
filenames <- list.files(pattern="\\.traces", full.names=TRUE)
gutlist <- lapply(filenames, read.neuron)
gutlistall <- as.neuronlist(gutlist)

Genotype1 = "DcrTripControl" #change for appropriate genotypes
Genotype2 = "DcrTraRNAi" #change for appropriate genotypes

#Run if subsetting gut by landmarks - otherwise skip to use full gut
#Generate landmarks file containing xyz coordinates of landmarks
gut <- read_csv(file=paste("Landmarks.csv",sep=""))
name <- gut[,1]
Landmark1 <- cbind(name,gut[,4:6])
Landmark2 <- cbind(name,gut[,8:10])
#Check landmarks and centrelines are in the same order
filenames == name
cbind(filenames,name)
#Subset centrelines by landmarks
gutlistsubset = neuronlist()
for(i in 1:dim(Landmark1)[1]) {
  distance <- rdist(as.matrix(Landmark1[i,2:4]),xyzmatrix(gutlistall[i]))
  minpoint1 <- which.min(distance)
  distance <- rdist(as.matrix(Landmark2[i,2:4]),xyzmatrix(gutlistall[i]))
  minpoint2 <- which.min(distance)
  gutlistsubset[[i]] = subset(gutlist[[i]],c(minpoint1:minpoint2))
}
gutlistall = gutlistsubset
#Plots to check subset by landmarks has worked
plot3d(gutlistsubset[[2]])
plot3d(gutlist[[2]], col='red')

#Make a csv file with the neuron names as the first column with the heading "File" and then subsequent columns as different the different info you want to record about the flies eg. sex, age
gutdataall <- read.csv("Trace_data_combined.csv", row.names="File") #load the csv file
#Add names to the gutlist
names(gutlistall) <- row.names(gutdataall)
# Add the csv as the dataframe for gutlist
data.frame(gutlistall) <- gutdataall
#CORRECT FOR UNIT DIFFERENCES BETWEEN DATASETS - converts all units to millimeters
UnitCorrection <- as.vector(gutdataall$Site)
UnitCorrection[UnitCorrection=='Site1'] <- 0.00296164
UnitCorrection[UnitCorrection=='Site2'] <- 0.002951345
UnitCorrection[UnitCorrection=='Site3'] <- 10
UnitCorrection[UnitCorrection=='Site4'] <- 25.4
UnitCorrection[UnitCorrection=='Site5'] <- 0.001
UnitCorrection <- as.numeric(UnitCorrection)
for (i in 1:length(UnitCorrection)){
  xyzmatrix(gutlistall[[i]]) <- xyzmatrix(gutlistall[[i]])*UnitCorrection[i] 
}

len <- vector()
gutlistre <- neuronlist()
for (i in 1:length(gutlistall)) {len[i]<-spine(gutlistall[[i]], rval = 'length')} # Finds all the spine lengths of every gut and puts them in a vector 'len'
for (i in 1:length(gutlistall)) {gutlistre[[i]] <- nat::resample(gutlistall[[i]], stepsize = len[i]/1000)}
gutlist1000 <- array(dim = c(2000,3,length(gutlistall)))
for (i in 1:length(gutlistall)) {gutlist1000[1:dim(xyzmatrix(gutlistre[[i]]))[1],,i] <- xyzmatrix(gutlistre[[i]])}
gutLandmark <- gutlist1000[1:1000,,]
#Using the gutLandmark file compares guts by 1000 pseudolandmarks

#INVERT data along x axis from bruker microCT
gutLandmarkinvert = gutLandmark
gutLandmarkinvert[,1,gutdataall$Site=='Site1'|gutdataall$Site=='Site2'|gutdataall$Site=='Site3'] = gutLandmark[,1,gutdataall$Site=='Site1'|gutdataall$Site=='Site2'|gutdataall$Site=='Site3']*-1
gutLandmark = gutLandmarkinvert

#Save gut length
gutLength <- as.data.frame(len)
row.names(gutLength) <- row.names(gutdataall)
write.csv(gutLength, "Gut_length.csv")

#check that filenames match between filename list and in csv trace data file
cbind(filenames,rownames(gutdataall)) 

#Extract data from gutdataall
ind <- row.names(gutdataall)
pops <- as.factor(gutdataall[,3])
sex <- as.factor(gutdataall[,1])
popssex <- as.factor(paste(gutdataall[,3],gutdataall[,1]))
batch <- as.factor(gutdataall[,18])

#Procrustes alignment
all.gpa <- gpagen(gutLandmark, PrinAxes = TRUE)

#PCA analysis
allPCA <-gm.prcomp(A=all.gpa$coords)
summary(allPCA)

#Generate images to visualise PCA axes
#Calculate consensus shape
consensus <- mshape(all.gpa$coords)

clear3d()
mfrow3d(1, 2, sharedMouse = TRUE)
#Plot hypothetical ends of PC1
plotRefToTarget(consensus,allPCA$shapes$shapes.PC1$min, method="points", mag=1)
plotRefToTarget(consensus,allPCA$shapes$shapes.PC1$max, method="points", mag=1)
#Save the plots as png and html
rgl.snapshot(filename = "PC1_minmax_all_usingPCAshapes.png")
writeWebGL(filename = "PC1_minmax_all_usingPCAshapes.html")

#Do the same for PC2
plotRefToTarget(consensus,allPCA$shapes$shapes.PC2$min, method="points", mag=1)
plotRefToTarget(consensus,allPCA$shapes$shapes.PC2$max, method="points", mag=1)
rgl.snapshot(filename = "PC2_minmax_all_usingPCAshapes.png")
writeWebGL(filename = "PC2_minmax_all_usingPCAshapes.html")

#Can also do the same for subsequent PCs by changing the number

#Plotting average shapes
consensustestV <- mshape(all.gpa$coords[,,gutdataall$Sex=='Virgin' & gutdataall$Genotype==Genotype1])
consensustestM <- mshape(all.gpa$coords[,,gutdataall$Sex=='Male' & gutdataall$Genotype==Genotype1])
consensusconV <- mshape(all.gpa$coords[,,gutdataall$Sex=='Virgin' & gutdataall$Genotype==Genotype2])
consensusconM <- mshape(all.gpa$coords[,,gutdataall$Sex=='Male' & gutdataall$Genotype==Genotype2])

clear3d()
mfrow3d(2, 2, sharedMouse = TRUE)
plot3d(consensusconM, axes = FALSE, xlab = "", ylab = "", zlab = "", xlim = NULL)
box3d()
aspect3d("iso")
plot3d(consensustestM, axes = FALSE, xlab = "", ylab = "", zlab = "", xlim = NULL)
box3d()
aspect3d("iso")
plot3d(consensusconV, axes = FALSE, xlab = "", ylab = "", zlab = "", xlim = NULL)
box3d()
aspect3d("iso")
plot3d(consensustestV, axes = FALSE, xlab = "", ylab = "", zlab = "", xlim = NULL)
box3d()
aspect3d("iso")
rgl.snapshot(filename = "AverageShapes_con_test_male_female.png")
rgl.postscript( filename= "AverageShapes_con_test_male_female.eps", fmt="eps", drawText=TRUE )
writeWebGL(filename = "AverageShapes_con_test_male_female.html")

#PCA plots
#Generate colours for PCA plot
col.pops <- rainbow(length(levels(popssex)))
names(col.pops) <- levels(popssex)
col.pops <- col.pops[match(popssex, names(col.pops))]

#make axes labels for PCA plot
xlab <- paste("Principal Component 2 ", "(", round(allPCA$pc.summary$importance[2,1]*100, 1), "%)", sep="")
ylab <- paste("Principal Component 3 ", "(", round(allPCA$pc.summary$importance[2,2]*100, 1), "%)", sep="")

GPA2d<-data.frame(allPCA$x)
PCA<-data.frame(GPA2d, popssex) #add sex/ other categories if present
fortify(PCA)
df<-PCA[c(1:120)] # (1:x) depends on how many PCs there are in the PCA

autoplot(prcomp(df), data = PCA, x=1, y=2, colour = "popssex", frame = TRUE, frame.type = 'norm', frame.colour = 'popssex', frame.alpha=0) + 
  scale_y_continuous(breaks = seq(-1, 1, by = 0.1), limits=c(-0.3,0.4)) + 
  scale_x_continuous(breaks = seq(-1, 1, by = 0.1), limits=c(-0.4,0.3)) + 
  theme(panel.background = element_rect(fill = NA, color = "black", size=1.5, linetype="solid"), panel.grid.major = element_line(size = 0.5, linetype = 'dashed',colour = "grey"), panel.grid.minor = element_blank(), axis.text=element_text(size=16), axis.title=element_text(size=18)) + 
  scale_color_discrete(name="Introgression Lines") + guides(fill = FALSE)

#Analysis
#The resultant Procrustes coordinates (found in all.gpa) can then be used for further statistical analysis.
#One example is Procrustes ANOVA with permutation procedures which quantifies the relative amount of shape variation attributable to one or more factors using a linear model and assesses this via permutation.
#Y is the shape data and x are one or more independent variable which can be discrete (as factors) or continuous.
#For these analyses, the data from GPA must be in a 2D array which can be done by geomorph.dataframe or two.d.array to convert the all.gpa$coords

gdf <- geomorph.data.frame(shape = all.gpa$coords, pops = pops, sex=sex, batch=batch)

#Analysis - Procrustes ANOVA with interaction term, with and without batch, typeIII
fit.pops <- procD.lm(shape ~ pops * sex + batch * pops + batch * sex, data = gdf, iter = 999, RRPP = TRUE, SS.type = "III")
fit.reduced <- procD.lm(shape ~ batch, data = gdf, iter = 999, RRPP = TRUE, SS.type = "III")
summary(fit.pops)
pwpops = pairwise(fit = fit.pops, fit.null = fit.reduced, groups = interaction(sex,pops))
summary(pwpops, test.type = "dist", confidence = 0.95, covariate = NULL)
summary(pwpops, test.type = "var", confidence = 0.95) #compares variation between the groups
fit.pops <- procD.lm(shape ~ pops * sex, data = gdf, iter = 999, RRPP = TRUE, SS.type = "III")
fit.reduced <- procD.lm(shape ~ 1, data = gdf, iter = 999, RRPP = TRUE, SS.type = "III")
summary(fit.pops)
pwpops = pairwise(fit = fit.pops, fit.null = fit.reduced, groups = interaction(sex,pops))
summary(pwpops, test.type = "dist", confidence = 0.95, covariate = NULL)
summary(pwpops, test.type = "var", confidence = 0.95) #compares variation between the groups
