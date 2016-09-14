# Hopefully we shall required not more than three seeds
# May be one seed during feature selection and 
# another seed during validation0
seed1 <- 100
seed2 <- 200
seed3 <- 300
source('Sepsis_MC_analysis_functions.R')

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
	if (is.null(alldata$dname) == FALSE) {
		dname <- alldata$dname
		if (dname == "AcbMero") {
			#Spart <- which(lclass %in% "S")
			#Rpart <- which(lclass %in% "R")
    		#SpartLen <- length(Spart)
    		#RpartLen <- length(Rpart)
			

			#Sodd <- seq(1, SpartLen, 2)
			#Seven <- seq(2, SpartLen, 2)
			#Rodd <- seq(1, RpartLen, 2)
			#Reven <- seq(2, RpartLen, 2)

			#SoddPart <- Spart[Sodd]
			#SevenPart <- Spart[Seven]
			#RoddPart <- Rpart[Rodd]
			#RevenPart <- Rpart[Reven]
			#trainPart <- c(SoddPart, RoddPart)
			#testpart <- c(SevenPart, RevenPart)
		
			#trainC <- lclass[trainPart]
			#testC <- lclass[testpart]

			#testC <- c(SevenPart, RevenPart)

			lclass <- alldata$lclass
			lstrain <- "RB197"
			lclass1 <- lclass[!(names(lclass) %in% lstrain)]
			odd <- seq(1, length(lclass1), 2)
    		even <- seq(2, length(lclass1), 2)
    		trainC <- lclass1[odd]
   			testC <- lclass1[even]
			
			testSpart <- subset(testC, testC == "S")
			testSpart.names <- names(testSpart)
			testSpart.vals <- as.character(testSpart)


			testRpart <- subset(testC, testC == "R")
            testRpart.names <- names(testRpart)
            testRpart.vals <- as.character(testRpart)

			testnames <- c(testSpart.names, lstrain, testRpart.names)
			testvals <- c(testSpart.vals, as.character(lclass[lstrain]), 
					testRpart.vals)
			names(testvals) <- testnames
			testC <- factor(testvals)


			alldata2 <- c(alldata, list(trainC = trainC, testC = testC))
			return (alldata2)

		}
	}
	
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


doFeatureSelection <- function(alldata, fsMethod, pval = 0.01) {
	if (fsMethod == "ReliefF") {
		alldata2 <- getFeaturesReliefF(alldata)
		return (alldata2)
	} else if (fsMethod == "rfRFE") {
		alldata2 <- getFeaturesRfRFE(alldata)
		return (alldata2)
	} else if (fsMethod == "pGreedy") {
		alldata2 <- getFeaturesPGreedy(alldata)
		return (alldata2)
	} else if (fsMethod == "SVMRFE") {
		alldata2 <- getFeaturesSVMRFE(alldata)
		return (alldata2)
	}
	else if (fsMethod == "Boruta") {
		alldata2 <- getFeaturesBoruta(alldata, pval) 
		return (alldata2)
	}
}


# We shall use all the features whatever returned by Boruta
getFeaturesBoruta <- function(alldata, pval) {
	trainC <- alldata$trainC
    ydata <- trainC
    cdata <- alldata$cdata
    xdata1 <- cdata[, names(trainC)]
    xdata2 <- t(xdata1)
	xdata <- data.frame(xdata2, trainC)

	BorutaRes <- Boruta(trainC ~., data = xdata, pValue = pval)
	finalDecision <- BorutaRes$finalDecision
	features <- names(finalDecision[finalDecision == "Confirmed"])
	alldata2 <- c(alldata, list(BorutaRes = BorutaRes, features = features))
    return (alldata2) 


}

getFeaturesSVMRFE <- function(alldata, ltimes = 5, featureCount = 2) {
	
	# create the index for five fold data with five times repeat
    trainC <- alldata$trainC
    ydata <- trainC
    cdata <- alldata$cdata
    xdata1 <- cdata[, names(trainC)]
    xdata <- t(xdata1)

    set.seed(seed1)
    finalCVCount <- getCVCount(trainC)
    index <- createMultiFolds(ydata, k = finalCVCount, times = ltimes)

	svmFuncs <- caretFuncs
	#svmFuncs$summary <- fivestats

	ctrl <- rfeControl(method = "repeatedcv",
		saveDetails = TRUE,
		number = finalCVCount,
		repeats = 5,
		returnResamp = "all",
		verbose = TRUE,
		functions = svmFuncs,
		index = index)


    varSeq <- seq(1, dim(xdata)[2] -1, by = 1)

	svmRFE <- rfe(x = xdata,
		y = ydata,
		sizes = varSeq,
		#metric = "ROC",
		rfeControl = ctrl,
		## Now options to train()
		method = "svmRadial",
		tuneLength = 12,
		preProc = c("center", "scale"),
		## Below specifies the inner resampling process
		trControl = trainControl(method = "cv",
		verboseIter = FALSE,
		classProbs = TRUE))

    features <- svmRFE$optVariables[1:featureCount]
    alldata2 <- c(alldata, list(svmRFE = svmRFE, features = features))
    return (alldata2) 
	
}

getCVCount <- function(classLabels) {

	standardCVCount <- 5
	smallClassCount <- min(table(classLabels))
	finalCVCount <- min(standardCVCount, smallClassCount)
	return (finalCVCount)
}


# this starts with seedNum probes with lowest p-values. The assumption is that there would be 
# one or two magic genes.

wilcoxon_single <- function(probe, res_part, sus_part) {
	res_single <- as.numeric(res_part[probe, ])
	sus_single <- as.numeric(sus_part[probe, ])

	lres <- wilcox.test(res_single, sus_single)
	pval <- lres$p.value
	return (pval)
}


ttest_single <- function(probe, res_part, sus_part) {
    res_single <- as.numeric(res_part[probe, ])
    sus_single <- as.numeric(sus_part[probe, ])

    lres <- t.test(res_single, sus_single)
    pval <- lres$p.value
    return (pval)
}



getPvalsAll <- function(ldata, lclass) {

    # Identify the pos genes and neg genes

    res_part <- ldata[, lclass == "R"]
    sus_part <- ldata[, lclass == "S"]

    probes <- rownames(res_part)

    pval_lst <- lapply(probes, wilcoxon_single, res_part, sus_part)
    #pval_lst <- lapply(probes, ttest_single, res_part, sus_part)
    pvals1 <- unlist(pval_lst)
    names(pvals1) <- probes
    pvals <- pvals1

}


drawPvals <- function (alldata, ltitle, plotname, outpath) {
	trainData <- data.frame(alldata$cdata[, names(alldata$trainC)])
    trainClass <- alldata$trainC

    trainpvals <- getPvalsAll(trainData, trainClass)

    testData <- data.frame(alldata$cdata[, names(alldata$testC)])
    testClass <- alldata$testC

    testpvals <- getPvalsAll(testData, testClass)

    train_pval_log2 <- log2(trainpvals)
    test_pval_log2 <- log2(testpvals)

	ldata <- data.frame(train_pval_log2, test_pval_log2)
	p <- ggplot(ldata, aes(train_pval_log2, test_pval_log2, label = names(trainpvals)))
	p + geom_point() + geom_text(hjust = 0, nudge_x = 0.05, size = 3) + ggtitle(ltitle)

	lfile <- paste0('Pvals_', ltitle)
	pdffile <- paste0(outpath, '/', lfile, '.pdf')
	ggsave(pdffile)


	ldata <- alldata$cdata
	lclass <- alldata$lclass

	allpvals <- getPvalsAll(ldata, lclass)
	all_pval_log2 <- log2(allpvals)

	ldata <- data.frame(train_pval_log2, all_pval_log2)
	p <- ggplot(ldata, aes(train_pval_log2, all_pval_log2, label = names(trainpvals)))
    p + geom_point() + geom_text(hjust = 0, nudge_x = 0.05, size = 3) + ggtitle(ltitle)

    lfile <- paste0('Pvals_all_', ltitle)
    pdffile <- paste0(outpath, '/', lfile, '.pdf')
    ggsave(pdffile)


}


getFeaturesPGreedy <- function(alldata, ltimes = 5, featureCount = 5, seedNum = 1) {

	trainData <- data.frame(alldata$cdata[, names(alldata$trainC)])
    trainClass <- alldata$trainC

    # Identify the pos genes and neg genes

    res_part <- trainData[, trainClass == "R"]
    sus_part <- trainData[, trainClass == "S"]

	probes <- rownames(res_part)
	
	pval_lst <- lapply(probes, wilcoxon_single, res_part, sus_part)
	#pval_lst <- lapply(probes, ttest_single, res_part, sus_part)
	pvals1 <- unlist(pval_lst)
	names(pvals1) <- probes
	pvals <- sort(pvals1)

	lseed <- pvals[1:seedNum]
	seed.genes <- names(lseed)

	# Now we have to apply another algorithm that would increase the set of features
	# without hurting the performance
	
	infoLst <- getInfo(alldata)
    disc_genes <- infoLst$disc_genes
    pos.genes <- infoLst$pos.genes
    neg.genes <- infoLst$neg.genes

	alldata2 <- getFeaturesGreedy(alldata, seed.genes, pos.genes, neg.genes,
        disc_genes)	
		
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

getFeaturesRfRFE <- function (alldata, ltimes = 5) {

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
		functions = newRF, index = index, verbose = TRUE) 
		  
	varSeq <- seq(10, dim(xdata)[2] -1, by = 2)
	
	#browser()
	rfRFE <- rfe(x = xdata, y = ydata, sizes = varSeq, metric = "ROC",
		rfeControl = ctrl, ntree = 1000)
	#features <- rfRFE$optVariables[1:featureCount]
	features <- rfRFE$optVariables
	alldata2 <- c(alldata, list(rfRFE = rfRFE, features = features))
	return (alldata2) 
}

# This would validate the data. Common for all the datasets, probably would 
# use random forest

validation <- function (alldata, lmethod, featureCount = 5) {

	testC <-  alldata$testC

	set.seed(seed2)
	finalCVCount <- getCVCount(testC)
	indexT <- createMultiFolds(testC, k = finalCVCount, times = 5)

	ctrlT <- trainControl(method = "repeatedcv", number = finalCVCount, 
		repeats = 5, returnResamp = "all", savePredictions = "all", 
		classProbs = TRUE, index = indexT)	

	features <- alldata$features[1:featureCount]
	cdata <- alldata$cdata
	xdata1 <- cdata[, names(testC)]
    xdata2 <- t(xdata1)
	testData <- data.frame(xdata2[, features], testC)
	modT <- train( testC ~ ., data = testData, trControl = ctrlT, method = lmethod)
	lresults <- modT$results	
	accuracy <- 0
	if (is.null(modT$bestTune$mtry) == FALSE) {
		bestMtry <- modT$bestTune$mtry
		accuracy <- lresults[lresults$mtry == bestMtry, "Accuracy"]	
	} else {
		accuracy <- lresults[, "Accuracy"]
	}
	alldata2 <- c(alldata, list(modT = modT, accuracy = accuracy))
	return (alldata2)

}

# Validation using full data; use train for training the model
# and test for prediction

validation_F <- function (alldata, lmethod, featureCount = 5) {

	set.seed(seed2)

	trainC <- alldata$trainC
	finalCVCount <- getCVCount(trainC)
    indexT <- createMultiFolds(trainC, k = finalCVCount, times = 5)

    ctrlT <- trainControl(method = "repeatedcv", number = finalCVCount,
        repeats = 5, returnResamp = "all", savePredictions = "all",
        classProbs = TRUE, index = indexT)


	features <- alldata$features[1:featureCount]
	cdata <- alldata$cdata
	xdata1 <- cdata[, names(trainC)]
    xdata2 <- t(xdata1)
	trainData <- data.frame(xdata2[, features], trainC)
	modT <- train( trainC ~ ., data = trainData, method = lmethod, 
			trControl = ctrlT)

	testC <-  alldata$testC
	ydata1 <- cdata[, names(testC)]
	ydata2 <- t(ydata1)
	testData <- data.frame(ydata2[, features])
	resClasses <- predict(modT, newdata = testData)
	resProbs <- predict(modT, newdata = testData, type = "prob")
	lmat <- confusionMatrix(resClasses, testC)
	accuracy <- as.numeric(lmat$overall["Accuracy"])

	alldata2 <- c(alldata, list(modT = modT, lmat = lmat, accuracy = accuracy, 
				resClasses = resClasses, resProbs = resProbs))
	
	return (alldata2)

}


plotCombinedMetric_F <- function(alldata, plotname, ltitle, outpath, 
						lmethod = "rf", featureCount = 5) {

	alldata2 <- validation_F(alldata, lmethod, featureCount)
    accuracy <- alldata2$accuracy
    lmetric <- getMetric_F(alldata2)


	resProbs <- alldata2$resProbs
	labels <- alldata$testC[rownames(resProbs)]
	predSample <- names(alldata$testC)
	predMIC <- as.numeric(alldata$MIC[predSample])
	lgroups1 <- paste0(predSample, "_", predMIC)
	lgroups <- factor(lgroups1, levels = lgroups1)

	mid_point <- max(alldata$MIC[alldata$lclass == 'S'])
    facet_var1 <- ifelse (predMIC <=mid_point , c('Sus'), c('Res'))
	facet_var <- factor(facet_var1, levels = c("Sus", "Res"))

	fMap <- alldata$fMap
    features <- alldata$features[1:featureCount]
    fVals1 <- fMap[features]
    fVals2 <- substr(fVals1, 1, 50)
    fVals <- paste(as.character(fVals2), collapse = "\n")

	ltitle1 <- paste0("Confidence of resistance\n", ltitle, ", 
		Accuracy = ", accuracy, "\nDecision score = ", lmetric, "\n", 
		fVals)
	
	probRes <- data.frame(lnames = rownames(resProbs), probs = resProbs[,1], 
			types = as.character(labels), lgroups = lgroups, 
			facet_var = facet_var)
	Palette1 <- c('red','forestgreen')
	plt <- ggplot(probRes, aes(x = lgroups, y = probs, colour = facet_var)) +
        geom_point(size = 3) +
		facet_grid(. ~ facet_var, scales = "free", space = "free") + 
		scale_colour_manual(values=Palette1) + 
		theme(axis.text.x = element_text(size=10,angle= 45))  + 
        xlab("Strain_MIC") + ylab("Probablity of resistance") +
		ggtitle(ltitle1) +  labs(colour='groups') 

    outpath <- lmethod
    lplotname <- paste0(ltitle,  ".pdf")
    lpoltpath = paste0(outpath, "/", lplotname)
    ggsave(lpoltpath)
	

}

plotCombinedMetric <- function(alldata, plotname, ltitle, outpath, featureCount = 5) {

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
	predLst <- lapply(lrnames, 
		function(x, pred2) {
			lobs <- as.character(pred2[x,"obs"])
			pred2[x, lobs]
		}, pred3
	)
	
	predObs <- unlist(predLst)
	predErr <- 1 - predObs
	pred4 <- data.frame(pred3, predErr)
	pred5 <- NULL
	if (is.null(modT$bestTune$mtry) == FALSE) {
		lmtry <- modT$bestTune$mtry
    	pred5 <- pred4[pred4$mtry == lmtry, ]
	} else {
		pred5 <- pred4
	}
    lgroups <- paste0(pred5$predSample, "_", pred5$predMIC)
    pred6 <- data.frame(pred5, lgroups)
    MICTest <- MIC[names(alldata$testC)]
    llevels <- paste0(names(MICTest), "_", MICTest)
    pred6$lgroups <- factor(pred6$lgroups, levels = llevels)

	# We take the logarithm
	# Distance compared to the middle point 2 in exponential
	
	exDistToMid <- pred6$predMIC/MICMid
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
    features <- alldata$features[1:featureCount]
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
	lpoltpath = paste0(outpath, "/", lplotname)
    ggsave(lpoltpath)

    #ltitle1 <- paste0("Prediction error\n", ltitle, ", Accuracy = ", accuracy, "\nDecision score = ", finalMetricStr, "\n", fVals)
	#plt <- ggplot(pred7, aes(x = lgroups, y = predErr)) + geom_boxplot(aes(fill=lcols[lgroups])) + 
	#	facet_grid(. ~ facet_var, scales = "free", space = "free") +  
	#	scale_colour_manual(values=Palette1) + xlab("Sample_MIC") + ylab("1 - Calling") +  
	#	ylim(0, 1) + theme(axis.text.x = element_text(size=10,angle= 45)) + 
	#	ggtitle(ltitle1)
	#lplotname <- paste0("Err_", plotname)
	#lpoltpath = paste0(outpath, "/", lplotname)
    #ggsave(lpoltpath)
    #ggsave(lplotname)
	return (plt_r)

}


getMetric_F <- function(alldata) {
	
	testC <- alldata$testC
	predSample <- alldata$testC
	resProbs <- alldata$resProbs

	lrnames <- rownames(resProbs)

	predLst <- lapply(lrnames, 
		function(x, resProbs, testC) {
			lclass1 <- as.character(testC[x])
			predVal <- resProbs[x, lclass1]
		}, resProbs, testC
	)
			
	predObs <- unlist(predLst)
	predErr <- 1 - predObs

	testMIC <- as.numeric(alldata$MIC[names(testC)])
	MICMid <- alldata$MICMid
	exDistToMid <- testMIC / MICMid
	
	logDistToMid <- log2(exDistToMid)
    corrDist <- logDistToMid
    corrDist[logDistToMid <= 0] <- abs(logDistToMid[logDistToMid <= 0]) + 1

	finalMetricAll <- predErr %*% corrDist
    lsum <- sum(corrDist)
    finalMetricAvg <- finalMetricAll / lsum
	return (finalMetricAvg)	
	
}



getMetric <- function(alldata) {
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
    predLst <- lapply(lrnames,
        function(x, pred2) {
            lobs <- as.character(pred2[x,"obs"])
            pred2[x, lobs]
        }, pred3
    )

    predObs <- unlist(predLst)
    predErr <- 1 - predObs
    pred4 <- data.frame(pred3, predErr)

	pred5 <- NULL
    if (is.null(modT$bestTune$mtry) == FALSE) {
        lmtry <- modT$bestTune$mtry
        pred5 <- pred4[pred4$mtry == lmtry, ]
    } else {
        pred5 <- pred4
    }

    lgroups <- paste0(pred5$predSample, "_", pred5$predMIC)
    pred6 <- data.frame(pred5, lgroups)
    MICTest <- MIC[names(alldata$testC)]
    llevels <- paste0(names(MICTest), "_", MICTest)
    pred6$lgroups <- factor(pred6$lgroups, levels = llevels)

    # We take the logarithm
    # Distance compared to the middle point 2 in exponential
   
    exDistToMid <- pred6$predMIC/MICMid
    logDistToMid <- log2(exDistToMid)
    corrDist <- logDistToMid
    corrDist[logDistToMid <= 0] <- abs(logDistToMid[logDistToMid <= 0]) + 1

    pred7 <- data.frame(pred6, corrDist)

    # Add one column for the facet_grid
    mid_point <- max(alldata$MIC[alldata$lclass == 'S'])
    facet_var <- ifelse (pred7$predMIC <= mid_point, c('Sus'), c('Res'))
    pred7$facet_var <- factor(facet_var, levels = c('Sus', 'Res'))

    finalMetricAll <- pred7$predErr %*% pred7$corrDist
    lsum <- sum(pred7$corrDist)
    finalMetricAvg <- finalMetricAll / lsum
   
    return (finalMetricAvg)

}

getFeaturesGreedy <- function(alldata, seed.genes, pos.genes, neg.genes, disc_genes, forwardThresh = 0.01, featureCount = 5) {

    yes.pos <- intersect(pos.genes, seed.genes)
    if (length(yes.pos) == 0) {
        yes.pos <- NULL
    }
    yes.neg <- intersect(neg.genes, seed.genes)
    if (length(yes.neg) == 0) {
        yes.neg = NULL
    }

    features <- forwardSearchWeighted_single(list(disc_genes), pos.genes, neg.genes, yes.pos=yes.pos, yes.neg=yes.neg, forwardThresh=forwardThresh, featureCount = featureCount)

    alldata2 <- c(alldata, list(features = unlist(features)))
}

getMetricSingle <- function(index, all.genes, pos.genes, neg.genes, 
		disc_genes, alldata2, lmethod) {
	parts <- strsplit(index, split = '_')
	parts_i <- unlist(parts)[1]
	parts_j <- unlist(parts)[2]
	i <- strtoi(parts_i)
	j <- strtoi(parts_j)
	seed.genes <- all.genes[c(i, j)]

    alldata3 <- getFeaturesGreedy(alldata2, seed.genes, pos.genes, neg.genes, disc_genes)
    alldata4 <- validation(alldata3, lmethod)
    accu <- alldata4$accuracy
    lMetric <- getMetric(alldata4)
	farr <- c(accu = accu, lMetric = lMetric, i = i, j = j)
	return(farr)

}

getInfo <- function(alldata) {
    disc_genes <- list()
    trainData <- data.frame(alldata$cdata[, names(alldata$trainC)])
    trainClass <- alldata$trainC
    disc_genes$genes <- data.frame(trainData)
    disc_genes$class <- as.integer(trainClass) -1

    # Identify the pos genes and neg genes

    res_part <- trainData[, trainClass == "R"]
    sus_part <- trainData[, trainClass == "S"]
    res_mean <- rowMeans(as.matrix(res_part))
    sus_mean <- rowMeans(as.matrix(sus_part))
    effect <- sus_mean - res_mean
    pos.genes <- names(effect[effect >= 0])
    neg.genes <- names(effect[effect < 0])
    retlst = list(disc_genes = disc_genes, pos.genes = pos.genes, neg.genes = neg.genes)
    return (retlst)
}

drawPlot <- function(index, all.genes, ldf, 
		alldata2, pos.genes, neg.genes, disc_genes, 
		dataname, partitionMethod, featureSelectionMethod, outpath, lmethod) {
	i <- ldf[index, "i"]
	j <- ldf[index, "j"]
	seed.genes <- all.genes[c(i, j)]

	alldata3 <- getFeaturesGreedy(alldata2, seed.genes, pos.genes, neg.genes,
    	disc_genes)
    # Do the validation and print the accuracy
    alldata4 <- validation(alldata3, lmethod)
    accu <- alldata4$accuracy
    lMetric <- getMetric(alldata4)
    # Generate a plot
    plotname <- paste0(dataname, "_", partitionMethod, "_", featureSelectionMethod, "_", i, "_", j, ".pdf")
    ltitle <- paste0(dataname, ", ", partitionMethod, ", ", featureSelectionMethod)
    plotCombinedMetric(alldata4, plotname, ltitle, outpath)
}

mainFuncGreedy <- function(datadir, dataFile, partitionMethod, lmethod = "rf", accuCutoff = 0.9, topPlotNum = 10) {
    load(dataFile)
    alldata2 <- doPartition(alldata, partitionMethod)

    all.genes <- rownames(alldata2$cdata)
    len.all.genes <- length(all.genes)
    infoLst <- getInfo(alldata2)
    disc_genes <- infoLst$disc_genes
    pos.genes <- infoLst$pos.genes
    neg.genes <- infoLst$neg.genes
    count <- 1
	featureSelectionMethod = "greedy"
	dataname <- strsplit(dataFile, split = "\\.")[[1]][1]
	outpath = paste0(dataname, "_", partitionMethod, "_", featureSelectionMethod)
	dir.create(file.path(datadir, outpath), showWarnings = FALSE)
	
    indices = c()
	for (i in c(1:len.all.genes)) {
        j = 1
        while (j < i) {
			part = paste0(i, '_', j)
			indices = c(part, indices)
			j = j+1
		}
	}

	res <- lapply(indices, getMetricSingle, all.genes, 
		pos.genes, neg.genes, disc_genes, alldata2, lmethod)

	ldf1 <- data.frame(do.call(rbind, res))
	ldf2 <- ldf1[order(-ldf1$accu, ldf1$lMetric),]

	# Now create plot for top ten
	ldf3 <- ldf2[1:topPlotNum,]
	lapply(1:topPlotNum, drawPlot, all.genes, ldf3, 
		alldata2, pos.genes, neg.genes, disc_genes,
		dataname, partitionMethod, featureSelectionMethod, outpath, lmethod)
}

getInfoFeatureWise <- function(featureCount, alldata, lmethod) {
	alldata2 <- validation(alldata, lmethod, featureCount)
	accuracy <- alldata2$accuracy
    lmetric <- getMetric(alldata2)

	return (c(accuracy = accuracy, lmetric = lmetric))

}


getInfoFeatureWise_F <- function(featureCount, alldata, lmethod) {
    alldata2 <- validation_F(alldata, lmethod, featureCount)
    accuracy <- alldata2$accuracy
    lmetric <- getMetric_F(alldata2)

    return (c(accuracy = accuracy, lmetric = lmetric))

}



drawProbPlotSpecific <- function(dataFile, partitionMethod, featureSelectionMethod, featureCount = 5, lmethod = "rf") {

	load(dataFile)
    alldata2 <- doPartition(alldata, partitionMethod)
    alldata3 <- doFeatureSelection(alldata2, featureSelectionMethod,
                pval = pval)
	dataname <- strsplit(dataFile, split = "\\.")[[1]][1]
    ltitle <- paste0("FeatureCount_", featureCount, "_", dataname, "_", 
				partitionMethod, "_", featureSelectionMethod)
	outpath <- lmethod
    lplotname <- paste0(ltitle,  ".pdf")
    lpoltpath = paste0(outpath, "/", lplotname)

	plotCombinedMetric_F(alldata3, lplotname, ltitle, outpath, lmethod, featureCount)

}

mainFuncFeatureWise_F <- function(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf") {

	load(dataFile)
    alldata2 <- doPartition(alldata, partitionMethod)
    alldata3 <- doFeatureSelection(alldata2, featureSelectionMethod, 
				pval = pval)

    featureLen <- length(alldata3$features)
    res1 <- lapply(2:featureLen, getInfoFeatureWise_F, alldata3, lmethod)
    res2 <- data.frame(featureCount = 2:featureLen, do.call(rbind, res1))
    names(res2)[2] <- "Accuracy"
    names(res2)[3] <- "Decision metric"
    res3 <- melt(res2, id.vars = 'featureCount')

    dataname <- strsplit(dataFile, split = "\\.")[[1]][1]
    ltitle <- paste0("Full_AllFeatures_", dataname, "_", partitionMethod, "_", featureSelectionMethod)
	 plt <- ggplot(res3, aes(x = featureCount, y = value, colour = variable)) +
        geom_point() + scale_x_continuous(breaks = seq(1, featureLen, 2)) +
        scale_y_continuous(breaks = seq(0, 1, 0.05)) +
        xlab("Number of features") + ylab("Accuracy/Decision metric") +
        ggtitle(ltitle)

    outpath <- lmethod
    lplotname <- paste0(ltitle,  ".pdf")
    lpoltpath = paste0(outpath, "/", lplotname)
    ggsave(lpoltpath)


}


mainFuncFeatureWise <- function(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf") {
	load(dataFile)
    alldata2 <- doPartition(alldata, partitionMethod)
    alldata3 <- doFeatureSelection(alldata2, featureSelectionMethod, pval = pval)

	featureLen <- length(alldata3$features)
	res1 <- lapply(1:featureLen, getInfoFeatureWise, alldata3, lmethod)
	res2 <- data.frame(featureCount = 1:featureLen, do.call(rbind, res1))
	names(res2)[2] <- "Accuracy"
	names(res2)[3] <- "Decision metric"
	res3 <- melt(res2, id.vars = 'featureCount')

	dataname <- strsplit(dataFile, split = "\\.")[[1]][1]
	ltitle <- paste0("AllFeatures_", dataname, "_", partitionMethod, "_", featureSelectionMethod)

	plt <- ggplot(res3, aes(x = featureCount, y = value, colour = variable)) + 
		geom_point() + scale_x_continuous(breaks = seq(1, featureLen, 2)) +
		scale_y_continuous(breaks = seq(0, 1, 0.05)) +  
		xlab("Number of features") + ylab("Accuracy/Decision metric") +
		ggtitle(ltitle)

    outpath <- lmethod
    lplotname <- paste0(ltitle,  ".pdf")
	lpoltpath = paste0(outpath, "/", lplotname)
    ggsave(lpoltpath)
	
}

mainFunc <- function(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf", pval = 0.01) {
	load(dataFile)
	alldata2 <- doPartition(alldata, partitionMethod)
	alldata3 <- doFeatureSelection(alldata2, featureSelectionMethod, pval = pval)
	#alldata4 <- validation_F(alldata3, lmethod, featureCount = 5)
	alldata4 <- validation(alldata3, lmethod)
	
	dataname <- strsplit(dataFile, split = "\\.")[[1]][1]
	plotname <- paste0(dataname, "_", partitionMethod, "_", featureSelectionMethod, ".pdf")
	ltitle <- paste0(dataname, ", ", partitionMethod, ", ", featureSelectionMethod)
	outpath <- lmethod
	dir.create(file.path(datadir, outpath), showWarnings = FALSE)

	plotCombinedMetric(alldata4, plotname, ltitle, outpath)

}


mainPvalPlotFunc <- function(dataFile, partitionMethod, lmethod = "rf") {
	load(dataFile)
    alldata2 <- doPartition(alldata, partitionMethod)
	dataname <- strsplit(dataFile, split = "\\.")[[r]][1]
    plotname <- paste0(dataname, "_", partitionMethod, ".pdf")
    ltitle <- paste0(dataname, ", ", partitionMethod)
    outpath <- lmethod
    dir.create(file.path(datadir, outpath), showWarnings = FALSE)
	drawPvals(alldata2, ltitle, plotname, outpath)

}	
