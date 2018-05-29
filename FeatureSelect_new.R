# Hopefully we shall required not more than three seeds
# May be one seed during feature selection and 
# another seed during validation0
seed1 <- 100
seed2 <- 200
seed3 <- 300
source('/home/nirmalya/research/featureselect/Sepsis_MC_analysis_functions.R')

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
	} 
}


getCVCount <- function(classLabels) {

	standardCVCount <- 5
	smallClassCount <- min(table(classLabels))
	finalCVCount <- min(standardCVCount, smallClassCount)
	return (finalCVCount)
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



# validation using the data on the feature selection space, using five fold
# cross validation
# The parameter datacov decides on the amount of data this is used 
# for validation

validation_FCV <- function(alldata, lmethod, featureCount = 5, datacov = 'testC') {

	testC <- NULL
	if (datacov == 'trainC') {
		testC <-  alldata$trainC
	} else if (datacov == 'testC') {
		testC <- alldata$testC
	} else if (datacov == 'allC') {
		testC <- alldata$lclass
	}

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


getInfoFeatureWise_FCV <- function(featureCount, alldata, lmethod) {
    alldata2 <- validation_FCV(alldata, lmethod, featureCount)
    accuracy <- alldata2$accuracy
    lmetric <- getMetric(alldata2)

    return (c(accuracy = accuracy, lmetric = lmetric))

}


mainFuncFeatureWise_FCV <- function(dataFile, partitionMethod, featureSelectionMethod, lmethod = "rf") {

	load(dataFile)
    alldata2 <- doPartition(alldata, partitionMethod)
    alldata3 <- doFeatureSelection(alldata2, featureSelectionMethod, 
				pval = pval)

    featureLen <- length(alldata3$features)
    res1 <- lapply(2:featureLen, getInfoFeatureWise_FCV, alldata3, lmethod)
    res2 <- data.frame(featureCount = 2:featureLen, do.call(rbind, res1))
    names(res2)[2] <- "Accuracy"
    names(res2)[3] <- "Decision metric"
    res3 <- melt(res2, id.vars = 'featureCount')

    dataname <- strsplit(dataFile, split = "\\.")[[1]][1]
    ltitle <- paste0("AllFeatures_FS_Test_", dataname, "_", partitionMethod, "_", featureSelectionMethod)
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


getDecisionPos <- function(R_probs, labels, C) {

    model <- list()
    labels_len <- length(labels)
    Q_dim <- 2 + labels_len
    Q <- matrix(0, Q_dim, Q_dim)
    Q[1,1] <- 1

    slack_mat <- diag(labels_len)
    ic_col <- rep(1, labels_len)
    A_part <- cbind(R_probs, ic_col)
    mult_var <- rep(0, labels_len)
    mult_var[labels == 'S'] <- -1
    mult_var[labels == 'R'] <- 1
    A_part2 <- mult_var * A_part
    A <- cbind(A_part2, slack_mat)

    rhs <- rep(1, labels_len)
    sense <- rep(">=", labels_len)

    lb <- c(rep(-Inf, 2), rep(0, labels_len))

    obj <- c(rep(0,2), rep(C, labels_len))

    model$Q <- Q
    model$A <- A
    model$obj <- obj
    model$lb <- lb
    model$sense <- sense
    model$rhs <- rhs

    res <- gurobi(model)
    final_val <- (-res$x[2]/res$x[1])
    return (final_val)

}

plotCombinedMetric_F <- function(alldata, plotname, ltitle, outpath,
                        lmethod = "rf", featureCount = 5, C = 100) {

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


    probRes <- data.frame(lnames = rownames(resProbs), probs = resProbs[,1],
            types = as.character(labels), lgroups = lgroups,
            facet_var = facet_var)
    Palette1 <- c('red','forestgreen')

    # Lets get the line

    R_probs <- resProbs[,1]
    final_pos <- getDecisionPos(R_probs, labels, C) 

    ltitle1 <- paste0("Confidence of resistance\n", ltitle, ",
        Accuracy = ", accuracy, "\nDecision score = ", lmetric, 
        "\nDecision position = ", final_pos, "\n",
        fVals)
    plt <- ggplot(probRes, aes(x = lgroups, y = probs, colour = facet_var)) +
        geom_point(size = 3) +
        facet_grid(. ~ facet_var, scales = "free", space = "free") +
        scale_colour_manual(values=Palette1) +
        theme(axis.text.x = element_text(size=10,angle= 45))  +
        geom_hline(yintercept = final_pos) +
        #theme(axis.text.x = element_blank(), axis.title.x = element_blank())  +
        xlab("Strain_MIC") + ylab("Probablity of resistance") +
        #ylab("Probablity of resistance") +
        ggtitle(ltitle1) +  labs(colour='groups')

    outpath <- lmethod
    lplotname <- paste0(ltitle,  ".pdf")
    lpoltpath = paste0(outpath, "/", lplotname)
    ggsave(lpoltpath)



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

