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

#l_dataFile <- 'KpCip3_allRespNormFold_runs1-2-4.RData'
l_dataFile <- 'KpMero3_2_fold.RData'
#l_dataFile <- 'AcbMero.RData'

lmethod = "rf"
fnum = 10


#get_features_on_all_data(dataFile = l_dataFile, featureSelectionMethod = "rfRFE", featureCount = fnum)
get_features_on_all_data(dataFile = l_dataFile, featureSelectionMethod = "ReliefF", featureCount = fnum)

#drawProbPlotSpecific(dataFile = l_dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", featureCount = fnum, lmethod = "rf") 
#drawProbPlotSpecific(dataFile = l_dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", featureCount = fnum, lmethod = "rf") 
#drawProbPlotSpecific(dataFile = l_dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", featureCount = fnum, lmethod = "rf") 
#drawProbPlotSpecific(dataFile = l_dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", featureCount = fnum, lmethod = "rf") 

#mainFuncFeatureWise_F(dataFile = 'KpCipII.RData', partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod = "rf")
#mainFuncFeatureWise_F(dataFile = 'KpCipII.RData', partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod = "rf")
#mainFuncFeatureWise_F(dataFile = 'KpCipII.RData', partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod = "rf")
#mainFuncFeatureWise_F(dataFile = 'KpCipII.RData', partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod = "rf")


#---------------------------------------------------------------------------
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)


#partitionMethod <- 'alternate'
#featureSelectionMethod <- 'rfRFE'

#mainFuncFeatureWise_F(dataFile = 'AcbMero.RData', partitionMethod, featureSelectionMethod, lmethod = "rf")
#mainFuncFeatureWise_F(dataFile = 'KpMero.RData', partitionMethod, featureSelectionMethod, lmethod = "rf")
#mainFuncFeatureWise_F(dataFile = 'KpGent.RData', partitionMethod, featureSelectionMethod, lmethod = "rf")

#mainFuncFeatureWise(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf")
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "Boruta", lmethod = "rf")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "Boruta", lmethod = "rf")

#dataFile <- 'KpCip.RData'
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "Boruta", lmethod = "rf", pval = 0.005)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "Boruta", lmethod = "rf")


#dataFile <- 'KpGent.RData'
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "Boruta", lmethod = "rf")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "Boruta", lmethod = "rf")


#dataFile <- 'KpMero.RData'
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "Boruta", lmethod = "rf")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "Boruta", lmethod = "rf")



#lmethod = 'rf'

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "pGreedy", lmethod)

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)


#dataFile <- 'KpCip.RData'

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "pGreedy", lmethod)

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

#dataFile <- 'KpGent.RData'

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "pGreedy", lmethod)

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

#dataFile <- 'KpMero.RData'
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "pGreedy", lmethod)

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE", lmethod)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE", lmethod)
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF", lmethod)
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF", lmethod)

