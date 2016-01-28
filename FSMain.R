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

S_index <- which(suscFac %in% 'S')
R_index <- which(suscFac %in% 'R')

split_f <- function(x, parts) {
	split(x, rep(1:length(parts), parts))
}

S_parts <- split_f(S_index, c(2, 2, 2, 3, 3))
R_parts <- split_f(R_index, c(2, 2, 2, 2, 2))
#S_parts[[1]] <- c(1,2,12)
#S_parts[[5]] <- c(10,11)
lEst <- "ReliefFexpRank"
#lEst <- "ReliefFbestK"
#lEst <- "ReliefFequalK"

# Five fold cross validation
all_index <- 1:length(suscFac)
tpred <- list()
true_test_c <- NULL
all_best_genes <- c()

for (j in 1:5) {

	test_index <- c(S_parts[[j]], R_parts[[j]])
	train_index <- setdiff(all_index, test_index)
    lsuscFac <- suscFac[train_index]
    tsuscFac <- suscFac[test_index]
    if (j == 1) {
		true_test_c <- tsuscFac
	} else {
		true_test_c <- as.factor(c(as.character(true_test_c), as.character(tsuscFac)))
	}
	train_data <- mydata4[train_index, ]
	test_data <- mydata4[test_index, ]
	ldata <- data.frame(train_data, lsuscFac)
	# get the set of features over the training samples
	estReliefF <- attrEval(lsuscFac ~ ., ldata, estimator= lEst)
	sortedFeatures <- sort(estReliefF, decreasing = TRUE)

    best_five <- sortedFeatures[1:5]
	best_ten <- sortedFeatures[1:10]
	print(best_five)
	all_best_genes <- c(all_best_genes, names(best_five))

	l1 <- 1
	for (k in 3:5) {

        best_k <- names(sortedFeatures[1:k])
		train_k <- data.frame(train_data[, best_k], lsuscFac)
		test_k <- data.frame(test_data[, best_k], tsuscFac)

       	# Now train a model and test it
		modelRF <- CoreModel(lsuscFac ~ ., train_k, model="rf",
			selectionEstimator="MDL",minNodeWeightRF=5,
			rfNoTrees=100, maxThreads=1)
		pred <- predict(modelRF, test_k, type="both")

		if (j == 1) {
			tpred[[l1]] <- pred
		} else {
			tpred[[l1]]$class <- factor(c(as.character(tpred[[l1]]$class), as.character(pred$class)))
			tpred[[l1]]$probabilities <- rbind(tpred[[l1]]$probabilities, pred$probabilities) 
		}

        l1 <- l1 + 1
    }
}


#library(ROCR)

#for (k in 1:3) {
#	xval <- tpred[[k]]$probabilities
#	yval <- true_test_c 
#	lpred <- prediction()
#}


mEvals <- list()
for (k in 1:3) {
	mEval <- modelEval(modelRF, true_test_c, tpred[[k]]$class, tpred[[k]]$prob)
	mEvals[[k]] <- mEval

	print(mEvals[[k]]$accuracy)
}

all_best_genes2 <- sort(all_best_genes, decreasing = FALSE)
print(table(all_best_genes2))


#estReliefF <- attrEval(suscFac ~ ., mydata5, estimator=lEst, ReliefIterations=30)

#estReliefF_S <- sort(estReliefF, decreasing = TRUE)
# Remove the rows that are useless
