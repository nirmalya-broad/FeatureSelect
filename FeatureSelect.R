# Hopefully we shall required not more than three seeds
# May be one seed during feature selection and 
# another seed during validation0
seed1 <- 100
seed2 <- 200
seed3 <- 300

# Partition the data

doPartition <- function (alldata, type) {
	if (type == "extreme") { 
		result <- doPartitionExtreme(alldata)
		return (result)
	} else if (type == "alternate") {
		result <- doPartitionAlternate(alldata)
		return (result)
	}
}

# Use the odd sampels for training and even one for testing/
# validation.

doPartitionAlternate <- function(alldata) {
	lclass <- alldata$lclass
	odd <- seq(1, length(lclass), 2)
	even <- seq(2, length(lclass), 2)
	trainC <- lclass[odd]
	testC <- lclass[even]
	alldata2 <- c(alldata, list(trainC = trainC, testC = testC))
	return (alldata2)
}

# Use the leftmost and rightmost section for feature selection
# the middle region for testing and validation

doPartitionExtreme <- function(alldata) {
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
	trainC <- factor(trainC1)

	testNames <- setdiff(names(lclass), names(trainC))
	testC <- factor(lclass[testNames])	
	
	alldata2 <- c(alldata, list(trainC = trainC, testC = testC))
    return (alldata2)

}


doFeatureSelection <- function(alldata, fsMethod) {
	if (fsMethod == "ReliefF") {
		alldata2 <- getFeaturesReliefF(alldata)
		return (alldata2)
	} else if (fsMethod == "rfRFE") {
		alldata2 <- getFeaturesRfRFE(alldata)
		return (alldata2)
	}
}

getCVCount <- function(classLabels) {

	standardCVCount <- 5
	smallClassCount <- min(table(classLabels))
	finalCVCount <- min(standardCVCount, smallClassCount)
	return (finalCVCount)
}

getFeaturesReliefF <- function(alldata, ltimes = 5, featureCount = 5, lEst = "ReliefFexpRank") {

	trainC <- alldata$trainC
    cdata <- alldata$cdata

	set.seed(seed1)
	finalCVCount <- getCVCount(trainC)
    index <- createMultiFolds(trainC, k = finalCVCount, times = ltimes)

	all_list <- c()
	for (j in 1:length(index)) {
		train_index <- index[[j]]
    	ydata <- trainC[train_index]
		xdata1 <- cdata[, names(ydata)]
		xdata2 <- t(xdata1)
		tdata <- data.frame(xdata2, ydata)
		estReliefF <- attrEval(ydata ~ ., tdata, estimator= lEst)
		all_list <- c(all_list, estReliefF)
		
	}

	all_list_df <- data.frame(names(all_list), as.numeric(all_list))
    colnames(all_list_df) <- c("Variables", "Values")
	finalList1 <- ddply(all_list_df,  .(Variables), function(x) mean(x$Values))
	names(finalList1)[2] <- "Values"
	finalList <- finalList1[order(finalList1$Values, decreasing = TRUE), ]
	final_genes <- as.character(finalList$Variables)
	
	features <- final_genes[1:featureCount]
	alldata2 <- c(alldata, list(features = features))	

}

getFeaturesRfRFE <- function (alldata, ltimes = 5, featureCount = 5) {

	# create the index for five fold data with five times repeat
	trainC <- alldata$trainC
	ydata <- trainC
	cdata <- alldata$cdata
	xdata1 <- cdata[, names(trainC)]
	xdata <- t(xdata1)
	
	set.seed(seed1)
	finalCVCount <- getCVCount(trainC)
	index <- createMultiFolds(ydata, k = finalCVCount, times = ltimes)

	newRF <- rfFuncs
	# Rerank is set true for reranking
	ctrl <- rfeControl(method = "repeatedcv", saveDetails = TRUE, 
		number = finalCVCount, repeats = 5, returnResamp = "all",  
		functions = newRF, index = index) 
		  
	varSeq <- seq(5, dim(xdata)[2] -1, by = 2)
	
	#browser()
	rfRFE <- rfe(x = xdata, y = ydata, sizes = varSeq, imetric = "ROC",
		rfeControl = ctrl, ntree = 1000)
	features <- rfRFE$optVariables[1:featureCount]
	alldata2 <- c(alldata, list(rfRFE = rfRFE, features = features))
	return (alldata2) 
}

# This would validate the data. Common for all the datasets, probably would 
# use random forest

validation <- function (alldata) {

	testC <-  alldata$testC

	set.seed(seed2)
	finalCVCount <- getCVCount(testC)
	indexT <- createMultiFolds(testC, k = finalCVCount, times = 5)

	ctrlT <- trainControl(method = "repeatedcv", number = finalCVCount, 
		repeats = 5, returnResamp = "all", savePredictions = "all", 
		classProbs = TRUE, index = indexT)	

	features <- alldata$features
	cdata <- alldata$cdata
	xdata1 <- cdata[, names(testC)]
    xdata2 <- t(xdata1)
	testData <- data.frame(xdata2[, features], testC)
	modT <- train( testC ~ ., data = testData, trControl = ctrlT)
	bestMtry <- modT$bestTune$mtry
	lresults <- modT$results	
	accuracy <- lresults[lresults$mtry == bestMtry, "Accuracy"]	
	alldata2 <- c(alldata, list(modT = modT, accuracy = accuracy))
	return (alldata2)

}

plotCombinedMetric <- function(alldata, plotname, ltitle) {

	modT <- alldata$modT
	testSample <- rownames(modT$trainingData)
    predSample <-  testSample[modT$pred$rowIndex]
    MIC <- alldata$MIC
	MICMid <- alldata$MICMid
	#MICMid <- 2
    pred1 <- modT$pred
	pred2  <- data.frame(predSample, pred1)	
	predMIC <- MIC[predSample]
	pred3  <-  data.frame(pred2, predMIC)

	# Now get the error with respect to observed; 
	lrnames <- rownames(pred3)
	predLst <- lapply(lrnames, 
		function(x, pred2) {
			lobs <- as.character(pred2[x,"obs"])
			pred2[x, lobs]
		}, pred3
	)
	
	predObs <- unlist(predLst)
	predErr <- 1 - predObs
	pred4 <- data.frame(pred3, predErr)
	lmtry <- modT$bestTune$mtry
    pred5 <- pred4[pred4$mtry == lmtry, ]
    lgroups <- paste0(pred5$predSample, "_", pred5$predMIC)
    pred6 <- data.frame(pred5, lgroups)
    MICTest <- MIC[names(alldata$testC)]
    llevels <- paste0(names(MICTest), "_", MICTest)
    pred6$lgroups <- factor(pred6$lgroups, levels = llevels)

	# We take the logarithm
	# Distance compared to the middle point 2 in exponential
	
	exDistToMid <- MICMid/pred6$predMIC
	logDistToMid <- log2(exDistToMid)
	corrDist <- logDistToMid
	corrDist[logDistToMid <= 0] <- abs(logDistToMid[logDistToMid <= 0]) + 1

	pred7 <- data.frame(pred6, corrDist)

	# Add one column for the facet_grid
	mid_point <- max(alldata$MIC[alldata$lclass == 'S'])
	facet_var <- ifelse (pred7$predMIC <=mid_point , c('Sus'), c('Res'))
	pred7$facet_var <- factor(facet_var, levels = c('Sus', 'Res'))

	finalMetricAll <- pred7$predErr %*% pred7$corrDist
	lsum <- sum(pred7$corrDist)
	finalMetricAvg <- finalMetricAll / lsum
	finalMetricStr <- sprintf("%0.4f", finalMetricAvg)

	fMap <- alldata$fMap
    features <- alldata$features
    fVals1 <- fMap[features]
	fVals2 <- substr(fVals1, 1, 50) 
    fVals <- paste(as.character(fVals2), collapse = "\n")

	resultsObs <- ddply(pred7, .(lgroups),
        function(x) {
            tableInfo = table(x$obs)
            countS = tableInfo["S"]
            countR = tableInfo["R"]
            countTotal = length(x$obs)
            if (countS == countTotal) {
                "Sus"
            } else if (countR == countTotal) {
                "Res"
            } else {
                "Mixed"
            }
        }
    )
	
	lcols <- resultsObs$V1
	names(lcols) <- resultsObs$lgroups
	Palette1 <- c('red','green','blue')
	accuracy1 <- alldata$accuracy
	accuracy <- sprintf('%0.4f', accuracy1)

    ltitle1 <- paste0("Confidence of resistance\n", ltitle, ", Accuracy = ", accuracy, "\nDecision score = ", finalMetricStr, "\n", fVals)
	plt_r <- ggplot(pred7, aes(x = lgroups, y = R)) + geom_boxplot(aes(fill=lcols[lgroups])) +
		facet_grid(. ~ facet_var, scales = "free", space = "free") +
		scale_colour_manual(values=Palette1) + xlab("Sample_MIC") +
		ylab("Probability of resistance") +  ylim(0, 1) +
		theme(axis.text.x = element_text(size=10,angle= 45)) + ggtitle(ltitle1)
	lplotname = paste0("Res_", plotname)
    ggsave(lplotname)

    ltitle1 <- paste0("Prediction error\n", ltitle, ", Accuracy = ", accuracy, "\nDecision score = ", finalMetricStr, "\n", fVals)
	plt <- ggplot(pred7, aes(x = lgroups, y = predErr)) + geom_boxplot(aes(fill=lcols[lgroups])) + 
		facet_grid(. ~ facet_var, scales = "free", space = "free") +  
		scale_colour_manual(values=Palette1) + xlab("Sample_MIC") + ylab("1 - Calling") +  
		ylim(0, 1) + theme(axis.text.x = element_text(size=10,angle= 45)) + 
		ggtitle(ltitle1)
	lplotname <- paste0("Err_", plotname)
    ggsave(lplotname)
	return (plt_r)

}

plotResults <- function(alldata, plotname, ltitle) {
	modT <- alldata$modT
	testSample <- rownames(modT$trainingData)
	predSample <-  testSample[modT$pred$rowIndex]
	MIC <- alldata$MIC
	pred1 <- modT$pred
	pred2  <- data.frame(predSample, pred1)
	predMIC <- MIC[predSample]
 	pred3  <-  data.frame(pred2, predMIC)
	maxPred <- unlist(lapply(pred3$R, function(x) {max(x, 1-x)}))

	pred4 <- data.frame(pred3, maxPred)
 	lmtry <- modT$bestTune$mtry
 	pred5 <- pred4[pred4$mtry == lmtry, ]
	lgroups <- paste0(pred5$predSample, "_", pred5$predMIC)
	pred6 <- data.frame(pred5, lgroups)
	MICTest <- MIC[names(alldata$testC)]
	llevels <- paste0(names(MICTest), "_", MICTest)
	pred6$lgroups <- factor(pred6$lgroups, levels = llevels)

	fMap <- alldata$fMap
	features <- alldata$features
	fVals1 <- fMap[features]
	fVals <- paste(as.character(fVals1), collapse = "\n")
	
	results1 <- ddply(pred6, .(lgroups), 
		function(x) {
			tableInfo = table(x$pred)
			countS = tableInfo["S"]
			countR = tableInfo["R"]
			countTotal = length(x$pred)
			if (countS == countTotal) {
				"Sus"
			} else if (countR == countTotal) {
				"Res"
			} else {
				"Mixed"
			}
		}
	)
	
	lcols <- results1$V1
	names(lcols) <- results1$lgroups
	accuracy <- alldata$accuracy
	ltitle1 <- paste0(ltitle, ", Accuracy = ", accuracy, "\n", fVals)
	Palette1 <- c('red','green','blue')
plt <- ggplot(pred6, aes(x = lgroups, y = maxPred)) + geom_boxplot(aes(fill=lcols[lgroups])) + scale_colour_manual(values=Palette1) + xlab("Sample_MIC") + ylab("Probability of calling") +  theme(axis.text.x = element_text(size=10,angle= 45)) + ggtitle(ltitle1)
	ggsave(plotname)

}


mainFunc <- function(dataFile, partitionMethod, featureSelectionMethod) {
	load(dataFile)
	alldata2 <- doPartition(alldata, partitionMethod)
	alldata3 <- doFeatureSelection(alldata2, featureSelectionMethod)
	alldata4 <- validation(alldata3)
	
	dataname <- strsplit(dataFile, split = "\\.")[[1]][1]
	plotname <- paste0(dataname, "_", partitionMethod, "_", featureSelectionMethod, ".pdf")
	ltitle <- paste0(dataname, ", ", partitionMethod, ", ", featureSelectionMethod)
	#plotResults(alldata4, plotname, ltitle)	
	plotCombinedMetric(alldata4, plotname, ltitle)

}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  #library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#library(caret)
#library(ggplot2)
#library(CORElearn)
#library(plyr)
#library(rentrez)
#library(grid)
#library(gridExtra)

#datadir <- '/home/nirmalya/research/DataDx'
#setwd(datadir)

#dataFile <- 'AcbMero.RData'

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE")
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF")


#dataFile <- 'KpCip.RData'

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE")
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF")

#dataFile <- 'KpGent.RData'

#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "rfRFE")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "rfRFE")
#mainFunc(dataFile, partitionMethod = "alternate", featureSelectionMethod = "ReliefF")
#mainFunc(dataFile, partitionMethod = "extreme", featureSelectionMethod = "ReliefF")


