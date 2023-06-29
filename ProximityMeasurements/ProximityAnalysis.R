#Package installation:
#install.packages("stringr")
#install.packages("RANN")
#install.packages("nat")
#install.packages("ggplot2")
#install.packages("tidyverse")
#install.packages("fields")


library(stringr)
library(RANN)
library(nat)
library(ggplot2)
library(tidyverse)
library(fields)

# Set up folder containing a folder with the list of segmentations of organX
# and a second folder with the list of segmentations for organY
setwd("~/DemoData_Proximity")  #master folder 

organX = "gut" 
organY = "ovary2"
distanceA = 0.5 #in MM change as appropriate - max distance organs are separated by

file_list_organX <- list.files(path="./Gut", pattern="*.obj", full.names=FALSE) 
file_list_organY<- list.files(path="./Ovary2", pattern="*.obj", full.names=FALSE)
fly = sapply(str_split(file_list_organX,"_"),head,1)

count = 0 # Count will say how many flies have organs X and Y at less than distance A
trueflies = matrix(, nrow=length(fly),ncol=1) # This will be a list of the flies where that is the case

# For each fly:
for(i in 1:length(fly)) {
  flyi = fly[i]
  # Load organ X
  for(j in 1:length(file_list_organX)) {
    if(grepl(fly[i], file_list_organX[j], fixed=TRUE)){
      organXi <- read.table(paste0("./Gut/",file_list_organX[j]),sep="",skip = 11)
      colnames(organXi) <- c("point","X","Y","Z")
      organXi <- organXi[organXi$point == "v",,]
      organXi[,2:4] <- organXi[,2:4]*0.844} #x0.844 to correct pixel calibration
  }
  # Load organ Y
  for(k in 1:length(file_list_organY)) {
    if(grepl(fly[i], file_list_organY[k], fixed=TRUE)){
      organYi <- read.table(paste0("./Ovary2/",file_list_organY[k]),sep="",skip = 11)
      colnames(organYi) <- c("point","X","Y","Z")
      organYi <- organYi[organYi$point == "v",,]
      organYi[,2:4] <- organYi[,2:4]*0.844}
  }
  # Find distances for every point in organ X to closest point in organ Y
  distances <- nn2(organYi[,2:4],query = organXi[,2:4], k=1)
  # Find if organ X is less than distance A to organ Y for current fly
  if(any(as.data.frame(distances[2]) < distanceA)){
    count = count+1
    trueflies[i] <- flyi
  }
}
trueflies <- na.omit(trueflies)

trace_dir = "./Centrelines"
# read gut traces into R as neurons. 
filenames <- list.files(path = trace_dir, pattern="\\.swc", full.names=TRUE)
filenames <- filenames[grepl("V", filenames, fixed=TRUE)]
centrelinenames <- sapply(str_split(filenames,"/"),tail,1)
centrelinenames <- sapply(str_split(centrelinenames,".s"),head,1)

#Ensure centrelines and segmentations are accessed in same order
fly == centrelinenames #check these are the same
#cbind(fly,centrelinenames) 

gutlist <- lapply(filenames, read.neuron)
gutlistall <- as.neuronlist(gutlist)

# split centreline into 100 equally spaced points for averaging
len <- vector()
gutlistre <- neuronlist()
for (i in 1:length(gutlistall)) {len[i]<-spine(gutlistall[[i]], rval = 'length')} # Finds all the spine lengths of every gut and puts them in a vector 'len'
for (i in 1:length(gutlistall)) {gutlistre[[i]] <- nat::resample(gutlistall[[i]], stepsize = len[i]/100)}
gutlist100 <- array(dim = c(200,3,length(gutlistall)))
for (i in 1:length(gutlistall)) {gutlist100[1:dim(xyzmatrix(gutlistre[[i]]))[1],,i] <- xyzmatrix(gutlistre[[i]])}
gutLandmark <- gutlist100[1:100,,]

gutDists <- matrix(, nrow=100, ncol=(length(fly)))
# For each fly:
for(i in 1:length(fly)) {
  flyi = fly[i]
  # Load organ X
  for(j in 1:length(file_list_organX)) {
    if(grepl(fly[i], file_list_organX[j], fixed=TRUE)){
      organXi <- read.table(paste0("./Gut/",file_list_organX[j]),sep="",skip = 11)
      colnames(organXi) <- c("point","X","Y","Z")
      organXi <- organXi[organXi$point == "v",,]
      organXi[,2:4] <- organXi[,2:4]*0.844}
  }
  # Load organ Y
  for(k in 1:length(file_list_organY)) {
    if(grepl(fly[i], file_list_organY[k], fixed=TRUE)){
      organYi <- read.table(paste0("./Ovary2/",file_list_organY[k]),sep="",skip = 11)
      colnames(organYi) <- c("point","X","Y","Z")
      organYi <- organYi[organYi$point == "v",,]
      organYi[,2:4] <- organYi[,2:4]*0.844}
  }
  # Find nearest 20 point on gut to each centreline point, and find the minimum distance of one of these points
  centreline_to_gut <- nn2(organXi[,2:4], query = gutLandmark, k=20)
  centreline_to_gut <- centreline_to_gut[1]
  for (ii in 1:100) {
    centre_to_gut_coords <- organXi[unlist(as.data.frame(centreline_to_gut)[ii,]),2:4]
    centreline_dists <- nn2(organYi[,2:4], query = centre_to_gut_coords, k=1)
    gutDists[ii,i] <- min(unlist(centreline_dists[2]))
  }
}

averageCentrelineDist <- rowMeans(gutDists)
stdevCentrelineDist <- apply(gutDists,1,sd)
ggplot()+geom_line(aes(x=c(1:100),y=as.vector(averageCentrelineDist))) + 
  geom_ribbon(aes(x=c(1:100),ymax=averageCentrelineDist + stdevCentrelineDist,ymin=averageCentrelineDist - stdevCentrelineDist), alpha = 0.4, fill = "grey70", linetype = 0) +
  scale_y_continuous(limits=c(-0.05,1), n.breaks = 5, expand = c(0,0)) +
  scale_x_continuous(limits=c(0,100), expand = c(0,0)) +
  theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
        axis.line=element_line(colour='grey20'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Extension cut by landmarks
gut <- read_csv(file=paste("Landmarks.csv",sep=""))
name <- gut[,1]
Landmark1 <- cbind(name,gut[,4:6])
Landmark2 <- cbind(name,gut[,8:10])
Landmark1 <- Landmark1[grepl("V", Landmark1$File, fixed=TRUE),]
Landmark2 <- Landmark2[grepl("V", Landmark2$File, fixed=TRUE),]

# Check filenames are the same and in the same order
fly == Landmark1$File
#cbind(fly,name)

# Find where landmark is on 100 point centreline
minpoint1 <- matrix(, nrow=(length(name)), ncol=1)
minpoint2 <- matrix(, nrow=(length(name)), ncol=1)
for(i in 1:dim(Landmark1)[1]) {
  distance <- rdist(as.matrix(Landmark1[i,2:4]),gutLandmark[,,i])
  minpoint1[i] <- which.min(distance)
  distance <- rdist(as.matrix(Landmark2[i,2:4]),gutLandmark[,,i])
  minpoint2[i] <- which.min(distance)
}

# Cut distance centrelines by landmarks and check whether it is closer than distance A to organY
countLandmarks = 0 # Count will say how many flies have landmark-cut gut and Y at less than distance A
truefliesLandmarks = matrix(, nrow=length(fly),ncol=1) # This will be a list of the flies where that is the case
for (i in 1:length(fly)) {
  if(any(gutDists[minpoint1[i]:minpoint2[i],i] < distanceA)){
    countLandmarks = countLandmarks+1
    truefliesLandmarks[i] <- fly[i]
  }
}
truefliesLandmarks <- na.omit(truefliesLandmarks)
truefliesLandmarks


