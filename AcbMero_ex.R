library(rgl)
#dataFile <- '/home/nirmalya/research/DataDx/AcbMero.RData'
#dataFile <- '/home/nirmalya/research/DataDx/KpGent.RData'
#dataFile <- '/home/nirmalya/research/DataDx/KpCip.RData'
dataFile <- '/home/nirmalya/research/DataDx/KpMero.RData'

load(dataFile)
cdata <- alldata$cdata
cdataT <- t(cdata)
#cdataTL <- log(cdataT)
cdataTL <- cdataT
lclass <- alldata$lclass
lnames <- names(lclass)
cdataTL.pca <- prcomp(cdataTL, center = FALSE, scale. = FALSE) 
#cdataTL.pca <- prcomp(cdataTL, center = TRUE, scale. = FALSE) 

scores <- data.frame(lclass, cdataTL.pca$x[,1:3])
llabels <- paste(as.character(lclass), names(lclass), sep = '_')
#set.seed(100)
#cl <- kmeans(scores[,2:4], 2)
#scores$cluster <- as.factor(cl$cluster)
#plot3d(scores[,2:4], col = scores$cluster, size =10)
#text3d(scores[,2:4], text=llabels, adj = c(1,1))
#plot3d(scores[,2:4], col = as.integer(lclass), size =10)
#text3d(scores[,2:4], text=llabels, adj = c(1,1))

parts <- strsplit(dataFile, "/")
llen <- length(parts[[1]])
ltitle <- parts[[1]][llen]
#library(rgl)

colors <- as.character(lclass)
colors[lclass == "R"] <- "black"
colors[lclass == "S"] <- "red"

plot3d(scores$PC1, scores$PC2, scores$PC3, col=colors, size =10, xlab = "PC1", ylab = "PC2", zlab = "PC3")
text3d(scores$PC1, scores$PC2, scores$PC3, text=llabels, adj = c(1,1))
title3d(ltitle)

plot(cdataTL.pca, type = "l")

summary(cdataTL.pca)


