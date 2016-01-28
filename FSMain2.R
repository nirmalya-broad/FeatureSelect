library(CORElearn)
# Load the dataset

# inFileRaw contains data without normalization
inFileRaw = '/home/nirmalya/research/DataDx/KpCip_for_Nirmalya_unscaled.txt'

# inFileProcessed contains data with normalization by Roby

inFileProcessed <- '/home/nirmalya/research/DataDx/KpCip_for_Nirmalya_scaled.txt'
inFile <- inFileProcessed

mydata <- read.table(inFile, sep = '\t', row.names = 1, header = TRUE, stringsAsFactors = FALSE)

temp <- mydata[, 12]
mydata[,12] <- mydata[,1]
mydata[,1] <- temp

#strain <- mydata[1, ]
Susc_Res <- mydata[1, ]
MIC <- mydata[2, ] 

metaRows <- 1:2

mydata2 <- mydata[-metaRows, ]
mydata3 <- t(data.matrix(mydata2))
#mydata4 <- mydata3[, !colSums(mydata3 <= 0)]
mydata4 <- mydata3
mydata4[mydata4 <= 0] <- NA
suscFac <- factor(unlist(Susc_Res))
#mydata5 <- data.frame(suscFac, mydata4)
mydata5 <- data.frame(mydata4, suscFac)


library(caret)

set.seed(100)
index <- createMultiFolds(suscFac, k = 5, times = 5)

lEst <- "ReliefFexpRank"

all_index <- 1:length(suscFac)

weights <- list()


for (j in 1:length(index)) {

	train_index <- index[[j]]
	test_index <- setdiff(all_index, train_index) 
    
	lsuscFac <- suscFac[train_index]
    tsuscFac <- suscFac[test_index]
	train_data <- mydata4[train_index, ]
	test_data <- mydata4[test_index, ]
	ldata <- data.frame(train_data, lsuscFac)
	# get the set of features over the training samples
	estReliefF <- attrEval(lsuscFac ~ ., ldata, estimator= lEst)
	sortedFeatures <- sort(estReliefF, decreasing = TRUE)

	weights[[j]] <- as.numeric(estReliefF)

}

weightsm <- do.call(rbind, weights)
av_weight <- colMeans(weightsm)
names(av_weight) <- colnames(mydata4)
s_weight <- sort(av_weight, decreasing = TRUE)

# Now classify the data using those five fold repeated cross validations

sgenes <- names(s_weight)

for (k in 2:10) {

	# top genes for this iteration
	tgenes <- sgenes[1:k]
	ac_count <- 0
	totalAc <- 0
	
	for (j in 1:length(index)) {
		train_index <- index[[j]]
    	test_index <- setdiff(all_index, train_index)

    	lsuscFac <- suscFac[train_index]
    	tsuscFac <- suscFac[test_index]
    	train_k <- data.frame(mydata4[train_index, tgenes], lsuscFac)
    	test_k <- data.frame(mydata4[test_index, tgenes], tsuscFac)
		modelRF <- CoreModel(lsuscFac ~ ., train_k, model="rf",
            selectionEstimator="MDL",minNodeWeightRF=5,
            rfNoTrees=100, maxThreads=1)
        pred <- predict(modelRF, test_k, type="both")

		mEval <- modelEval(modelRF, tsuscFac, pred$class, pred$probabilities)
		totalAc <- totalAc + mEval$accuracy * length(test_index)
		ac_count <- ac_count + length(test_index)
    }
	av_accu <- totalAc / ac_count

	cat("Accuracy with top ", k, "genes: ", av_accu, "\n")
}


for (1:10)
cat("Top ten genes: ", sgenes[1:10], "\n")
