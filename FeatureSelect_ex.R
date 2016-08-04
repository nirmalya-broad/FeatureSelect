doFourPartition <- function(alldata) {

	lclass <- alldata$lclass
    Spart <- subset(lclass, lclass == "S")
    Rpart <- subset(lclass, lclass == "R")

    SpartLen <- length(Spart)
    RpartLen <- length(Rpart)

    leftTrain1 <- Spart[1:ceiling(SpartLen/2)]
    leftTrain <- as.character(leftTrain1)
    names(leftTrain) <- names(leftTrain1)
   
    rightTrain1 <- Rpart[ceiling(RpartLen/2 + 1):RpartLen]
    rightTrain <- as.character(rightTrain1)
    names(rightTrain) <- names(rightTrain1)
    trainC1 <- c(leftTrain, rightTrain)
    extremeClass <- factor(trainC1)

    testNames <- setdiff(names(lclass), names(extremeClass))
    interClass <- factor(lclass[testNames])

	# Now divide the train and test parts into two halves each.
	
	doPartsAlternate <- function(lclass) {
		odd <- seq(1, length(lclass), 2)
		even <- seq(2, length(lclass), 2)
		partA <- lclass[odd]
		partB  <- lclass[even]
		llist <- list(partA = partA, partB = partB)
		return (llist)
	}

	extremeParts <- doPartsAlternate(extremeClass)
	interParts <- doPartsAlternate(interClass)
	extremeA <- extremeParts$partA
	extremeB <- extremeParts$partB

	interA <- interParts$partA
	interB <- interParts$partB

	alldata2 <- c(alldata, list(extremeA = extremeA, extremeB = extremeB, 
			interA = interA, interB = interB))
	
}

prepareBoxData <- function(alldata, datatype, featureSize = 5) {

	modT <- alldata$modT
    testSample <- rownames(modT$trainingData)
    predSample <-  testSample[modT$pred$rowIndex]
    MIC <- alldata$MIC
    MICMid <- alldata$MICMid
    pred1 <- modT$pred
    pred2  <- data.frame(predSample, pred1)
    predMIC <- MIC[predSample]
    pred3  <-  data.frame(pred2, predMIC)

    # Now get the error with respect to observed;
    lrnames <- rownames(pred3)
    predLst <- lapply(lrnames, function(x, pred2) { 
		lobs <- as.character(pred2[x,"obs"])
		pred2[x, lobs]}, 
		pred3)
   
    predObs <- unlist(predLst)
    predErr <- 1 - predObs
    pred4 <- data.frame(pred3, predErr)
    #lmtry <- modT$bestTune$mtry
	lmtry <- featureSize
    pred5 <- pred4[pred4$mtry == lmtry, ]
    lgroups <- paste0(pred5$predSample, "_", pred5$predMIC)
    pred6 <- data.frame(pred5, lgroups)
    MICTest <- MIC[names(alldata$testC)]
    llevels <- paste0(names(MICTest), "_", MICTest)
    pred6$lgroups <- factor(pred6$lgroups, levels = llevels)

	# A good idea may be to add a column indicating whether the data is 
	# similar of different.

	pred7 <- data.frame(pred6, datatype)

	return (pred7)

}


doPart <- function(alldata, featurePart, similar, different, fsMethod, 
	plotname, ltitle) {

	# Now we shall extract the features from FSPart and validate on both 
	# similar and different. At present we shall validate on similar and 
	# different separately. One suggestion is to combine similar and 
	# different together. We shall do it, if it does not work the other
	# way.

	alldata$trainC <- featurePart
	alldata2 <- doFeatureSelection(alldata, fsMethod)
	
	# Part A -  validate on the same type of data
	alldataA <- alldata2
	alldataA$testC <- similar
	alldataA1 <- validation(alldataA)

	# Part B - validate on the other type of data
	alldataB <- alldata2
	othersC <- unlist(list(different[[1]], different[[2]]))
	othersC1 <- sort(othersC, decreasing = TRUE)
	alldataB$testC <- othersC1
	alldataB1 <- validation(alldataB)

	# Now we have to combine the results and print on a plot. 
	# On the boxplot we have to specifically show that there are 
	# two different kinds of data, same vs different (e.g. 
	# susceptible vs resistent) using two different colors.
	# Also on the plot there has to be a vertical line dividing
	# the two classes.

	boxDataSimilar <- prepareBoxData(alldataA1, datatype = "similar")
	boxDataDifferent <- prepareBoxData(alldataB1, datatype = "different")
	boxAll <- rbind(boxDataSimilar, boxDataDifferent)
	
	# Have to sort the data with respect to proper MIC. That would 
	# Help us to understand it clarly and also to organize similar and
	# different dataset.

	lgroups <- boxAll$lgroups
	llevels <- levels(lgroups)
	lMIC <- as.numeric(unlist(lapply(strsplit(llevels, "_"), 
		function(x) {
			x[[2]]
		}
	)))

	lMICdf <- data.frame(llevels, lMIC)
	lMICdf2 <- lMICdf[with(lMICdf, order(lMIC)),]
	llevels_sorted <- as.character(lMICdf2$llevels)	
	lgroups2 <- factor(lgroups, levels = llevels_sorted)
	boxAll2 <- boxAll
	boxAll2$lgroups <- lgroups2

	facet_var <- ifelse (boxAll2$predMIC <=2 , c('Sus'), c('Res'))
    boxAll2$facet_var <- factor(facet_var, levels = c('Sus', 'Res'))


	fMap <- alldata2$fMap
    features <- alldata2$features
    fVals1 <- fMap[features]
    fVals2 <- substr(fVals1, 1, 50)
    fVals <- paste(as.character(fVals2), collapse = "\n")


	ltitle1 <- paste0("Confidence on different regions\n", ltitle, "\n", fVals)
	
	plt_r <- ggplot(boxAll2, aes(x = lgroups, y = R)) + 
		geom_boxplot(aes(fill = datatype)) + 
		facet_grid(. ~ facet_var, scales = "free", space = "free") +
		xlab("Sample_MIC") + ylab("1 - Calling") +  
        ylim(0, 1) + theme(axis.text.x = element_text(size=10,angle= 45)) +
        ggtitle(ltitle1)
	
	lplotname <- paste0("Comp_", plotname)
	ggsave(lplotname)
	return (plt_r)

}

doCore <- function(alldata, FSPart, validationPart, fsMethod, 
	plotnamePart, fsType) {

	plotnamePart2 = paste(plotnamePart, fsMethod, fsType, sep = "_")
	titlePart = paste(plotnamePart, fsMethod, fsType, sep = ",")	
	
	FSPartA <- FSPart$partA
	FSPartB <- FSPart$partB

	fsPart = "A"
	plotname = paste0(plotnamePart2, "_", fsPart, ".pdf")
	ltitle = paste0(titlePart, ",", fsPart) 
	doPart(alldata, FSPartA, FSPartB, validationPart, fsMethod, plotname, 
		ltitle)

	fsPart = "B"
	plotname = paste0(plotnamePart2, "_", fsPart, ".pdf")
	ltitle = paste0(titlePart, ",", fsPart) 
	doPart(alldata, FSPartB, FSPartA, validationPart, fsMethod, plotname, 
		ltitle)

}

doPatterns <- function(dataFile, fsMethod) {

	load(dataFile)

	alldata2 <- doFourPartition(alldata)
	extremeA <- alldata2$extremeA
	extremeB <- alldata2$extremeB
	interA <- alldata2$interA
	interB <- alldata2$interB

	plotnamePart <- strsplit(dataFile, split = "\\.")[[1]][1]


	# Pattern I

	FSPart <- list(partA = extremeA, partB = extremeB)
	validationPart <- list(partA = interA, partB = interB)
	fsType = "extreme"
	doCore(alldata, FSPart, validationPart, fsMethod, plotnamePart, fsType)

	# Pattern II

	FSPart <- list(partA = interA, partB = interB)
	validationPart <- list(partA = extremeA, partB = extremeB)
	fsType = "intermediate"
	doCore(alldata, FSPart, validationPart, fsMethod, plotnamePart, fsType)

}


#library(caret)
#library(ggplot2)
#library(CORElearn)
#library(plyr)
#library(rentrez)
#library(grid)
#library(gridExtra)


#source('/home/nirmalya/research/featureselect/FeatureSelect.R')
#datadir <- '/home/nirmalya/research/DataDx'
#setwd(datadir)


#dataFile <- 'AcbMero.RData'
#doPatterns(dataFile, fsMethod = "rfRFE")
#doPatterns(dataFile, fsMethod = "ReliefF")

#dataFile <- 'KpCip.RData'
#doPatterns(dataFile, fsMethod = "rfRFE")
#doPatterns(dataFile, fsMethod = "ReliefF")

#dataFile <- 'KpGent.RData'
#doPatterns(dataFile, fsMethod = "rfRFE")
#doPatterns(dataFile, fsMethod = "ReliefF")
