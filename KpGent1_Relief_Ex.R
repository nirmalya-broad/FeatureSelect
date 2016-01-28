library(caret)
library(CORElearn)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(plyr)

datadir <- '/home/nirmalya/research/DataDx/'
outDir <- paste0(datadir, '/KpGent1_Relief_Ex/')
dir.create(outDir, showWarnings = FALSE)
datafile <- paste0(datadir, '/forNirmalya_KpGent1_norm.txt')
micfile <- paste0(datadir, '/forNirmalya_KpGent1_MICs.txt')

expData1 <- read.table(datafile, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
expData <- data.frame(t(expData1))
MICData <- read.table(micfile, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

MIC1 <-as.matrix(MICData)
MIC2 <- MIC1[order(MIC1),]

catego <- function(arr) {
	ltrainY1 <- arr 
	names(ltrainY1) <- paste0(names(ltrainY1), "_gent")
	ltrainY <- ltrainY1
	ltrainY[ltrainY1 <=2] = "S"
	ltrainY[ltrainY1 > 2] = "R"
	ltrainY <- factor(ltrainY)
	return (ltrainY)

}

part1 <- MIC2[1:11]
part2 <- MIC2[27:31]
ex_part <- c(part1, part2)
m_part <- MIC2[12:26]

trainY <- catego(ex_part)
testY <- catego(m_part)


# The data is now ready 

xdata1 <- expData
ydata <- trainY
xdata <- xdata1[names(trainY), ]

testX <- xdata1[names(testY), ]

######## Do some expr #####

index <- createMultiFolds(trainY, k = 5, times = 5)
fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...))
newRF <- rfFuncs
ctrl <- rfeControl(method = "repeatedcv", saveDetails = TRUE, number = 5, repeats = 5, returnResamp = "all",  verbose = TRUE, functions = newRF, index = index)
 # verbose = TRUE, functions = newRF, index = index)
varSeq <- seq(1, length(xdata)-1, by = 2)

rfRFE <- rfe(x = xdata, y = ydata, sizes = varSeq, metric = "ROC",
   rfeControl = ctrl,
    ## now pass options to randomForest()
    ntree = 1000)

indexT <- createMultiFolds(testY, k = 5, times = 5)

ctrlT <- trainControl(method = "repeatedcv", number = 5, repeats = 5,  verbose = TRUE, index = indexT)

testVars <- rfRFE$optVariables
testData <- data.frame(testX[, testVars], testY)
modT <- train(testY ~., data = testData, trControl = ctrlT)
modT
##---------Common part, copied from other file ------------------##

getFeatures <- function(xdata, ydata, index) {
	lEst <- "ReliefFexpRank"

	all_index <- 1:length(ydata)

	all_list <- c()
	for (j in 1:length(index)) {

    	train_index <- index[[j]]
    	test_index <- setdiff(all_index, train_index)

    	lydata <- ydata[train_index]
		cat("j: ", tydata, "\n")
    	tydata <- ydata[test_index]
    	train_data <- xdata[train_index, ]
    	test_data <- xdata[test_index, ]
    	ldata <- data.frame(train_data, lydata)
    	# get the set of features over the training samples
    	estReliefF <- attrEval(lydata ~ ., ldata, estimator= lEst)

    	all_list <- c(all_list, estReliefF)

	}

	all_list_df <- data.frame(names(all_list), as.numeric(all_list))
	colnames(all_list_df) <- c("Variables", "Values")
	return(all_list_df)

}


all_index <- 1:length(ydata)
index <- createMultiFolds(trainY, k = 5, times = 5)
all_list_df <- getFeatures(xdata, ydata, index)

#index <- createMultiFolds(suscFac, k = 5, times = 5)

finalList1 <- ddply(all_list_df,  .(Variables), function(x) mean(x$Values))
names(finalList1)[2] <- "Values" 
finalList <- finalList1[order(finalList1$Values, decreasing = TRUE), ]
final_genes <- as.character(finalList$Variables)


# Now classify the data using those five fold repeated cross validations

xdata <- xdata <- xdata1[names(testY), ]
ydata <- testY


all_index <-1:length(testY) 
index <- createMultiFolds(testY, k = 5, times = 5)

TestFunc <- function(xdata, ydata, final_genes, index) {
	R_pos <- which(ydata %in% 'R')
	S_pos <- which(ydata %in% 'S')

	lmat_lst <- list()
	setwd(outDir)
	for (k in 2:10) {

		# top genes for this iteration
		ac_count <- 0
		totalAc <- 0

		lmat <- matrix(0, length(index), length(all_index))
		rownames(lmat) <- names(index)
		colnames(lmat) <- names(ydata)

		for (j in 1:length(index)) {
			tgenes <- final_genes[1:k]
			train_index <- index[[j]]
			test_index <- setdiff(all_index, train_index)

			lydata <- ydata[train_index]
			tydata <- ydata[test_index]
			train_k <- data.frame(xdata[train_index, tgenes], lydata)
			test_k <- data.frame(xdata[test_index, tgenes], tydata)
			modelRF <- CoreModel(lydata ~ ., train_k, model="rf",
				selectionEstimator="MDL",minNodeWeightRF=5,
				rfNoTrees=100, maxThreads=1)
			pred <- predict(modelRF, test_k, type="both")

			mEval <- modelEval(modelRF, tydata, pred$class, pred$probabilities)
			totalAc <- totalAc + mEval$accuracy * length(test_index)
			if (mEval$accuracy < 1)
			  cat ("j pos: ", j , "\n")
			ac_count <- ac_count + length(test_index)

			# fill up the lmat for the jth training/test partition
			# 1: for training S
			# 2: for training R
			# 3: TS - true susceptible 
			# 4: TR - true resistant
			# 5: FS - false susceptible
			# 6: FR - false resistant
			# This all should serve the purpose.
			
			train_S <- intersect(S_pos, train_index)
			train_R <- intersect(R_pos, train_index)
			test_class <- ydata[test_index]
			pred_class <- pred$class
			TS <- test_index[which(test_class == 'S' & pred_class == 'S')]
			TR <- test_index[which(test_class == 'R' & pred_class == 'R')]
			FS <- test_index[which(test_class == 'R' & pred_class == 'S')]
			FR <- test_index[which(test_class == 'S' & pred_class == 'R')]

			lmat[j, train_S] <- "Train Sus"
			lmat[j, train_R]  <- "Train Res"
			lmat[j, TS] <-  "True Sus"
			lmat[j, TR] <- "True Res"
			lmat[j, FS] <- "False Sus"
			lmat[j, FR] <- "False Res"
		}
		av_accu <- totalAc / ac_count

		lmat_lst[[k]] <- lmat
		cat("Accuracy with top ", k, "genes: ", av_accu, "\n")
		
		#col<-c("grey","green", "white","blue","black", "red")
		#image(1:dim(lmat)[2],1:dim(lmat)[1], t(lmat), col= col)
		#lmat1 <- lmat[, MIC_Ordered]
		lmat1 <- lmat
		lmat2 <- melt(lmat1)

		myColors <- brewer.pal(length(levels(factor(lmat2[,3]))),"Dark2")
		names(myColors) <- levels(lmat2$values)
		colScale <- scale_fill_manual(name = "Group",values = myColors)

		file_name <- paste0("Plot_", k, ".pdf")
		#cols <- c("black", "gray", "blue", "green", "red", "yellow")
		ggplot(lmat2, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + theme(axis.text.x = element_text(size=8,angle=90)) + colScale
		ggsave(file_name)


	}

		return (lmat_lst)
}

llist <- TestFunc(xdata, ydata, final_genes, index)

