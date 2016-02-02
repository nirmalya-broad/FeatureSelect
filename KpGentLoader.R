# This file would convert the KpCip dataset into an R object.
# For each dataset, we may have a separate data loader. If 
# there are only few data format we can generalize them.

#args <- commandArgs(trailingOnly = TRUE)
#print (args)

#datadir <- args[1]
#infile <- args[2]
#outdir <- args[3]
 
datadir <- '/home/unix/nirmalya/Desktop/DataDx2'
infile <- 'forNirmalya_KpGent1_norm.txt'
MICFile <- 'forNirmalya_KpGent1_MICs.txt'
outfile <- 'KpGent.RData'
miccoundary <- 2

inpath <- paste0(datadir, '/', infile)
expdata1 <- read.csv(inpath, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)

MICPath <- paste0(datadir, '/', MICFile)
MICdata1 <- read.table(MICPath, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
MICdata2 <- MICdata1[,1]
names(MICdata2) <- rownames(MICdata1)
MIC <- MICdata2[order(MICdata2)]

cnames1 <- names(expdata1)
cnames <- gsub("_gent$", "", cnames1)
names(expdata1) <- cnames

expdata2 <- expdata1[, names(MIC)]
expdata3 <- lapply(expdata2, function(x) as.numeric(x))
cdata <- do.call('cbind', expdata3)
rownames(cdata) <- rownames(expdata2)

lclass1 <- rep("", length(MIC))
names(lclass1) <- names(MIC)
lclass1[MIC <= miccoundary] = "S"
lclass1[MIC > miccoundary] = "R"
lclass <- factor(lclass1)

alldata <- list(lclass = lclass, MIC = MIC, cdata = cdata)

outpath <- paste0(datadir, '/', outfile)
save(alldata, file = outfile)
