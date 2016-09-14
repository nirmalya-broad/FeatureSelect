library(caret)
library(ggplot2)
library(CORElearn)
library(plyr)
library(rentrez)
library(grid)
library(gridExtra)
library(FSelector)
library(Biocomb)
library(Boruta)
library(reshape2)

source('FeatureSelect.R')

datadir <- '/home/nirmalya/research/DataDx'
setwd(datadir)

dataFile <- 'KpMero.RData'

drawProbPlotSpecific(dataFile = "KpMero.RData", partitionMethod = "alternate", featureSelectionMethod = "rfRFE", featureCount = 5, lmethod = "rf") 


mainFuncFeatureWise_F(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf")
mainFuncFeatureWise(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf")
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "Boruta", lmethod = "rf")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "Boruta", lmethod = "rf")

dataFile <- 'KpCip.RData'
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "Boruta", lmethod = "rf", pval = 0.005)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "Boruta", lmethod = "rf")


dataFile <- 'KpGent.RData'
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "Boruta", lmethod = "rf")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "Boruta", lmethod = "rf")


dataFile <- 'KpMero.RData'
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "Boruta", lmethod = "rf")
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "Boruta", lmethod = "rf")



lmethod = 'rf'

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "pGreedy", lmethod)

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)


dataFile <- 'KpCip.RData'

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "pGreedy", lmethod)

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

dataFile <- 'KpGent.RData'

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "pGreedy", lmethod)

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

dataFile <- 'KpMero.RData'
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "pGreedy", lmethod)

mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

