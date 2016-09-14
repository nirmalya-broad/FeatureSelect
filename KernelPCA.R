datadir <- '/home/nirmalya/research/DataDx'
setwd(datadir)

library(kernlab)
library(ggplot2)

dataFile <- 'KpMero.RData'
load(dataFile)
ldata1 <- alldata$cdata
ldata2 <- data.frame(t(ldata1))

kpc <- kpca(~ ., data = ldata2, kernel = "rbfdot", kpar = list(sigma = lsigma),
            features = 3)   
kpc <- kpca(~ ., data = ldata2, kernel = "polydot", kpar = list(degree = 3),
kpc <- kpca(~ ., data = ldata2, kernel = "splinedot", kpar = list(), features = 3)   

lsigma <- 0.05
kpc <- kpca(~ ., data = ldata2, kernel = "rfdot", 
		kpar = list(sigma = lsigma), features = 3)   

comps1 <- rotated(kpc)
comps2 <- data.frame(comps1)
colnames(comps2) <- c("PC1", "PC2", "PC3")
lclass1 <- alldata$lclass
lclass1['RB270'] <- 'R'
lclass1['RB412'] <- 'S'
#lclass1['RB403'] <- 'R'

labels1 <- rownames(comps2)
MIC1 <- as.character(alldata$MIC)
labels2 <- paste0(labels1, '_', MIC1)

p <- ggplot(comps2, aes(PC1, PC2, label = labels2))
#p <- ggplot(comps2, aes(PC2, PC3, label = labels2))
p + geom_point(aes(colour = lclass1)) + 
	geom_text(hjust = 0, nudge_x = 0.05, size = 3) 
#ggsave('KpMero_RBF_0.05_PCA.pdf')

                                                                                  




