library(CORElearn)

# Load the dataset

# inFileRaw contains data without normalization
inFileRaw = '/home/nirmalya/research/DataDx/KpCip_for_Nirmalya_unscaled.txt'

# inFileProcessed contains data with normalization by Roby

inFileProcessed <- '/home/nirmalya/research/DataDx/KpCip_for_Nirmalya_scaled.txt'
inFile <- inFileProcessed

mydata <- read.table(inFile, sep = '\t', row.names = 1, header = TRUE, stringsAsFactors = FALSE)

#strain <- mydata[1, ]
Susc_Res <- mydata[1, ]
MIC <- mydata[2, ] 

metaRows <- 1:2

mydata2 <- mydata[-metaRows, ]
mydata3 <- t(data.matrix(mydata2))
mydata4 <- mydata3[, !colSums(mydata3 <= 0)]
suscFac <- factor(unlist(Susc_Res))
mydata5 <- data.frame(suscFac, mydata4)

estReliefF <- attrEval(suscFac ~ ., mydata5, estimator="ReliefFexpRank", ReliefIterations=30)

estReliefF_S <- sort(estReliefF, decreasing = TRUE)
# Remove the rows that are useless
