library(caret)
library(ggplot2)
library(CORElearn)
library(plyr)
library(rentrez)
library(grid)
library(gridExtra)


source('/home/nirmalya/research/featureselect/FeatureSelect.R')
source('/home/nirmalya/research/featureselect/FeatureSelect_ex.R')
datadir <- '/home/nirmalya/research/DataDx'
setwd(datadir)


dataFile <- 'AcbMero.RData'
doPatterns(dataFile, fsMethod = "rfRFE")
doPatterns(dataFile, fsMethod = "ReliefF")

dataFile <- 'KpCip.RData'
doPatterns(dataFile, fsMethod = "rfRFE")
doPatterns(dataFile, fsMethod = "ReliefF")

dataFile <- 'KpGent.RData'
doPatterns(dataFile, fsMethod = "rfRFE")
doPatterns(dataFile, fsMethod = "ReliefF")


