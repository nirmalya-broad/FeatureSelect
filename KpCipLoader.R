# This file would convert the KpCip dataset into an R object.
# For each dataset, we may have a separate data loader. If 
# there are only few data format we can generalize them.

#args <- commandArgs(trailingOnly = TRUE)
#print (args)

#datadir <- args[1]
#infile <- args[2]
#outdir <- args[3]
 
datadir <- '/home/unix/nirmalya/Desktop/DataDx2'
infile <- '/KpCip_for_Nirmalya_scaled.txt'
outfile <- 'KpCip.RData'

inpath <- paste0(datadir, '/', infile)
expdata1 <- read.csv(inpath, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

class1 <- expdata1[1,]
class2 <- as.character(class1)
names(class2) <- names(class1)
lclass <- factor(class2)

MIC1 <- expdata1[2, ]
MIC <- as.numeric(MIC1)
names(MIC) <- names(MIC1)

cdata1 <- expdata1[-c(1,2), ]
cdata2 <- lapply(cdata1, function(x) as.numeric(x))
cdata <- do.call (cbind, cdata2)
rownames(cdata) <- rownames(cdata1)

alldata <- list(lclass = lclass, MIC = MIC, cdata = cdata)

outpath <- paste0(datadir, '/', outfile)
save(alldata, file = outfile)
