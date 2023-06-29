#Input: files from FIJI macro 'Macro_measureIntensityAlong2DGutCentreline.ijm (provided in DemoData folder)
#Output: 'AllProfiles' table of position values normalised for gut length against intensity, images concatenated
#Output2: 'LMAdjProfiles' table of position values adjusted for 1 landmark against intensity, images concatenated

#Package installation:
#install.packages("tidyverse")
#install.packages("fields")

library(tidyverse)
library(fields)

#change working directory to folder where files from FIJI macro are saved
setwd("~/DemoData_HandtraRNAiBnl")

#Opens centreline files, adds columns for name and sex (sex must be written in filename!)
centrelinefiles <- list.files(path=getwd(), pattern="\\Centreline.csv", full.names=FALSE)
gutcentreline <- read_csv(file=centrelinefiles[1])
name <- strsplit(centrelinefiles[1], "_Intensity") %>% sapply('[',1)
if (grepl("VirginMale",name)){
  sex = rep(c("VirginMale"),dim(gutcentreline)[1])
}
if (grepl("VirginFemale", name)){
  sex = rep(c("VirginFemale"),dim(gutcentreline)[1])
}
if (grepl("DicertraTRiP_Pos",name)){
  geno = rep(c("dcrtraRNAi"),dim(gutcentreline)[1])
}
if (grepl("DicerCherryTRiP_Pos", name)){
  geno = rep(c("TripControl"),dim(gutcentreline)[1])
}
sexGenotype = paste0(sex, geno)

namelist <- rep(c(name),dim(gutcentreline)[1])
AllProfiles <- cbind(gutcentreline,sexGenotype,namelist)
LMAdjProfiles <- cbind(gutcentreline,sexGenotype,namelist)

cllength <- gutcentreline[dim(gutcentreline)[1],2]

#finds position of landmark and adjusts position values
xy <- read_csv(file=paste0(name,"_xyCoordinatesOfLengthandDistances.csv"))
Landmark <- read_csv(file=paste0(name,"_Landmarks.csv"))
distance <- rdist(as.matrix(Landmark[1,2:3]),xy[,2:3])
minpoint <- which.min(distance)
toLMlenxy <- sum(xy[1:minpoint,4])/1.2168  #divide by pixel calibration
Landmarkpos <- which.min(abs(gutcentreline$Position-toLMlenxy))
Landmarkdis <- gutcentreline$Position[Landmarkpos]
LMAdjProfiles$Position[1:Landmarkpos] <- as.numeric(LMAdjProfiles$Position[1:Landmarkpos]) / as.numeric(Landmarkdis)*100
LMAdjProfiles$Position[(Landmarkpos+1):length(LMAdjProfiles$Position)] <- (((as.numeric(LMAdjProfiles$Position[(Landmarkpos+1):length(LMAdjProfiles$Position)])-Landmarkdis) / as.numeric((cllength - Landmarkdis)))*100+100)

#normalises for length without landmark adjustment
AllProfiles$Position <- as.numeric(AllProfiles$Position) / as.numeric(cllength)*100

#repeats for rest of files in folder and concatenates
for (i in 2:length(centrelinefiles)){
  gutcentreline <- read_csv(file=centrelinefiles[i])
  name <- strsplit(centrelinefiles[i], "_Intensity") %>% sapply('[',1)
  if (grepl("VirginMale",name)){
    sex = rep(c("VirginMale"),dim(gutcentreline)[1])
  }
  if (grepl("VirginFemale", name)){
    sex = rep(c("VirginFemale"),dim(gutcentreline)[1])
  }
  if (grepl("DicertraTRiP_Pos",name)){
    geno = rep(c("dcrtraRNAi"),dim(gutcentreline)[1])
  }
  if (grepl("DicerCherryTRiP_Pos", name)){
    geno = rep(c("TripControl"),dim(gutcentreline)[1])
  }
  sexGenotype = paste0(sex, geno)
  
  namelist <- rep(c(name),dim(gutcentreline)[1])
  profile <- cbind(gutcentreline,sexGenotype,namelist)
  adjprofiles <- cbind(gutcentreline,sexGenotype,namelist)
  cllength <- gutcentreline[dim(gutcentreline)[1],2]
  
  #finds position of landmark and adjusts position values
  xy <- read_csv(file=paste0(name,"_xyCoordinatesOfLengthandDistances.csv"))
  Landmark <- read_csv(file=paste0(name,"_Landmarks.csv"))
  distance <- rdist(as.matrix(Landmark[1,2:3]),xy[,2:3])
  minpoint <- which.min(distance)
  toLMlenxy <- sum(xy[1:minpoint,4])/1.2168  #divide by pixel calibration
  Landmarkpos <- which.min(abs(gutcentreline$Position-toLMlenxy))
  Landmarkdis <- gutcentreline$Position[Landmarkpos]
  adjprofiles$Position[1:Landmarkpos] <- as.numeric(adjprofiles$Position[1:Landmarkpos]) / as.numeric(Landmarkdis)*100
  adjprofiles$Position[(Landmarkpos+1):length(adjprofiles$Position)] <- (((as.numeric(adjprofiles$Position[(Landmarkpos+1):length(adjprofiles$Position)])-Landmarkdis) / as.numeric((cllength - Landmarkdis)))*100+100)
  LMAdjProfiles <- rbind(LMAdjProfiles,adjprofiles)

  #normalises for length without landmark adjustment
  profile$Position <- as.numeric(profile$Position) / as.numeric(cllength)*100
  AllProfiles <- rbind(AllProfiles,profile)
}

write.csv(AllProfiles,file = "AllProfiles.csv")
write.csv(LMAdjProfiles,file = "LMAdjProfiles.csv")

#plot - mean and sd for landmark adjusted
binsize = 5
VirginControlPlot <- LMAdjProfiles %>% filter(sexGenotype == "VirginFemaleTripControl") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
MaleControlPlot <- LMAdjProfiles %>% filter(sexGenotype == "VirginMaleTripControl") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
VirginTraPlot <- LMAdjProfiles %>% filter(sexGenotype == "VirginFemaledcrtraRNAi") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
MaleTraPlot <- LMAdjProfiles %>% filter(sexGenotype == "VirginMaledcrtraRNAi") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
all <- bind_rows(list("VirginControl" = VirginControlPlot,"MaleControl" = MaleControlPlot,"VirginTra"=VirginTraPlot, "MaleTra" =MaleTraPlot),.id = "id")
eb <- aes(ymax = Intensity_mean + Intensity_sd, ymin = Intensity_mean - Intensity_sd)
ggplot(data = all, aes(x = Position_mean, y = Intensity_mean, color=id, fill = id)) + 
  geom_line(size = 2) + 
  geom_ribbon(eb, alpha = 0.4, fill = "grey70", linetype = 0) +
  scale_y_continuous(limits=c(-50,100), n.breaks = 5, expand = c(0,0)) +
  scale_x_continuous(limits=c(0,200), expand = c(0,0)) +
  theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
      panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
      axis.line=element_line(colour='grey20'),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position="none")

#new plot - mean and sd for non-landmark adjusted
binsize = 2
VirginControlPlot <- AllProfiles %>% filter(sexGenotype == "VirginFemaleTripControl") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
MaleControlPlot <- AllProfiles %>% filter(sexGenotype == "VirginMaleTripControl") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
VirginTraPlot <- AllProfiles %>% filter(sexGenotype == "VirginFemaledcrtraRNAi") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
MaleTraPlot <- AllProfiles %>% filter(sexGenotype == "VirginMaledcrtraRNAi") %>% subset(select=c(Position,Intensity)) %>% group_by(Position_bin = cut(Position, seq(0, 200, binsize),include.lowest = TRUE)) %>% summarise_all(list(mean=mean,sd=sd))
all <- bind_rows(list("VirginControl" = VirginControlPlot,"MaleControl" = MaleControlPlot,"VirginTra"=VirginTraPlot, "MaleTra" =MaleTraPlot),.id = "id")
eb <- aes(ymax = Intensity_mean + Intensity_sd, ymin = Intensity_mean - Intensity_sd)
ggplot(data = all, aes(x = Position_mean, y = Intensity_mean, color=id, fill = id)) + 
  geom_line(size = 2) + 
  geom_ribbon(eb, alpha = 0.4, fill = "grey70", linetype = 0) +
  scale_y_continuous(limits=c(-50,100), n.breaks = 5, expand = c(0,0)) +
  scale_x_continuous(limits=c(0,100), expand = c(0,0)) +
  theme(panel.background=element_rect(fill='white'), plot.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'), panel.grid.minor.y = element_blank(),
        axis.line=element_line(colour='grey20'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="none")
