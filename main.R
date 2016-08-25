library(caret)
library(ggplot2)
library(CORElearn)
library(plyr)
library(rentrez)
library(grid)
library(gridExtra)

source('FeatureSelect.R')

datadir <- '/home/nirmalya/research/DataDx'
setwd(datadir)

dataFile <- 'AcbMero.RData'
lmethod = "gaussprLinear"
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)


dataFile <- 'KpCip.RData'

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

dataFile <- 'KpGent.RData'

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

dataFile <- 'KpMero.RData'
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

