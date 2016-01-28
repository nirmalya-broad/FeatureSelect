library(caret)
library(CORElearn)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

datadir <- '/home/nirmalya/research/DataDx/'
outDir <- paste0(datadir, '/KpGent1_Relief/')
dir.create(outDir, showWarnings = FALSE)
datafile <- paste0(datadir, '/forNirmalya_KpGent1_norm.txt')
micfile <- paste0(datadir, '/forNirmalya_KpGent1_MICs.txt')

expData1 <- read.table(datafile, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

MICData <- read.table(micfile, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

MIC1 <-as.matrix(MICData)
MIC2 <- MIC1[order(MIC1),]
MIC_Ordered <- paste0(names(MIC2), "_gent")

expData <- data.frame(t(expData1))


suscNames <- rownames(subset(MICData, MIC <= 2))
suscIds <- paste0(suscNames, '_gent')
resNames <- rownames(subset(MICData, MIC > 2))
resIds <- paste0(resNames, '_gent')

yvalS <- rep("S", length(suscIds))
names(yvalS) <- suscIds

yvalR <- rep("R", length(resIds))
names(yvalR) <- resIds

yvals1 <- c(yvalS, yvalR)

yvals <- factor(yvals1)

# The data is now ready 

xdata1 <- expData
ydata <- yvals
xdata <- xdata1[names(ydata), ]


##---------Common part, copied from other file ------------------##


all_index <- 1:length(ydata)
index <- createMultiFolds(all_index, k = 5, times = 5)

#index <- createMultiFolds(suscFac, k = 5, times = 5)

suscFac <- ydata
mydata4 <- xdata

lEst <- "ReliefFexpRank"

all_index <- 1:length(suscFac)

weights <- list()

totalAc <- 0
ac_count <- 0

s_list <- list()
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

    weights[[j]] <- as.numeric(estReliefF)

	s_weight <- sort(estReliefF, decreasing = TRUE)
	sgenes <- names(s_weight)
	s_list[[j]] <- sgenes

}

# Now classify the data using those five fold repeated cross validations

R_pos <- which(ydata %in% 'R')
S_pos <- which(ydata %in% 'S')

setwd(outDir)
for (k in 2:10) {

    # top genes for this iteration
    ac_count <- 0
    totalAc <- 0

    lmat <- matrix(0, length(index), length(all_index))
    rownames(lmat) <- names(index)
    colnames(lmat) <- names(ydata)

    for (j in 1:length(index)) {
		sgenes <- s_list[[j]]
    	tgenes <- sgenes[1:k]
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

    cat("Accuracy with top ", k, "genes: ", av_accu, "\n")
	
    #col<-c("grey","green", "white","blue","black", "red")
	#image(1:dim(lmat)[2],1:dim(lmat)[1], t(lmat), col= col)
    lmat1 <- lmat[, MIC_Ordered]
	lmat2 <- melt(lmat1)

    myColors <- brewer.pal(length(levels(factor(lmat2[,3]))),"Dark2")
    names(myColors) <- levels(lmat2$values)
    colScale <- scale_fill_manual(name = "Group",values = myColors)

    file_name <- paste0("Plot_", k, ".pdf")
    #cols <- c("black", "gray", "blue", "green", "red", "yellow")
    ggplot(lmat2, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + theme(axis.text.x = element_text(size=8,angle=90)) + colScale
	ggsave(file_name)

}


