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

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE")
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF")


dataFile <- 'KpCip.RData'

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE")
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF")

dataFile <- 'KpGent.RData'

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE")
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF")

