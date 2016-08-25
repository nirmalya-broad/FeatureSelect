library(caret)
library(ggplot2)
library(CORElearn)
library(plyr)
library(rentrez)
library(grid)
library(gridExtra)
library(rmeta)
library(multtest)

source('FeatureSelect.R')
source('Sepsis_MC_analysis_functions.R')

datadir <- '/home/nirmalya/research/DataDx'
setwd(datadir)

dataFile <- 'AcbMero.RData'
#dataFile <- 'KpMero.RData'
#dataFile <- 'KpGent.RData'
#dataFile <- 'KpCip.RData'
load(dataFile)

mainFuncGreedy(datadir, dataFile, partitionMethod = "alternate", accuCutoff = 0.95)



