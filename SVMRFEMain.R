#library(CORElearn)
library(caret)

set.seed(100)
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
mydata4 <- mydata3[, !colSums(mydata3 <= 0)]
#mydata4 <- mydata3
#mydata4[mydata4 <= 0] <- NA
suscFac <- factor(unlist(Susc_Res))
#mydata5 <- data.frame(suscFac, mydata4)
#mydata5 <- data.frame(mydata4, suscFac)

xdata <- data.frame(mydata4)
ydata <- suscFac

fiveStats <- function(...) c(twoClassSummary(...), defaultSummary(...))
newRF <- rfFuncs
newRF$summary <- fiveStats
index <- createMultiFolds(ydata, k = 5, times = 5)

varSeq <- seq(1, length(xdata)-1, by = 2)

ctrl <- rfeControl(method = "repeatedcv", number = 5, repeats = 5, 
  functions = newRF, index = index)
 # verbose = TRUE, functions = newRF, index = index)

rfRFE <- rfe(x = xdata, y = ydata, sizes = varSeq, metric = "ROC", 
   rfeControl = ctrl, 
    ## now pass options to randomForest()
    ntree = 1000)

#set.seed(100)


# Print the 

ctrl <- rfeControl(method = "repeatedcv", number = 5, repeats = 5, functions = svmFuncs, index = index)

svmRFE <- rfe(x = xdata, y = ydata, sizes = varSeq,
  metric = "ROC", rfeControl = ctrl, ## Now options to train()
  method = "svmRadial", tuneLength = 12,
  preProc = c("center", "scale"),
  ## Below specifies the inner resampling process
  trControl = trainControl(method = "cv",
  verboseIter = FALSE, classProbs = TRUE))



