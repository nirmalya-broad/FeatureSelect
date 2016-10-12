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

source('FeatureSelect_new.R')

datadir <- '/home/nirmalya/research/DataDx'
setwd(datadir)

partitionMethod <- 'alternate'
featureSelectionMethod <- 'rfRFE'
lmethod <- 'rf'


dataFile <- 'KpMero.RData'
drawProbPlotSpecific(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", featureCount = 5, lmethod = "rf")


dataFile <- 'AcbMero.RData'
drawProbPlotSpecific(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", featureCount = 5, lmethod = "rf")


dataFile <- 'KpGent.RData'
drawProbPlotSpecific(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", featureCount = 5, lmethod = "rf")

dataFile <- 'KpMero.RData'
mainFuncFeatureWise_FCV(dataFile, partitionMethod, featureSelectionMethod, lmethod)

dataFile <- 'AcbMero.RData'
mainFuncFeatureWise_FCV(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf")

dataFile <- 'KpGent.RData'
mainFuncFeatureWise_FCV(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf")





