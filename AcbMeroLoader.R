# This file would convert the KpCip dataset into an R object.
# For each dataset, we may have a separate data loader. If 
# there are only few data format we can generalize them.

#args <- commandArgs(trailingOnly = TRUE)
#print (args)

#datadir <- args[1]
#infile <- args[2]
#outdir <- args[3]

library(rentrez) 
datadir <- '/home/nirmalya/research/DataDx'
infile <- 'AcbMero1_norm_foldInd.csv'
MICFile <- 'AcbMero1_strains_MICs.txt'
outfile <- 'AcbMero.RData'
miccoundary <- 2
MICMid <- 2

inpath <- paste0(datadir, '/', infile)
expdata1 <- read.csv(inpath, header = TRUE, row.names = 1, stringsAsFactors = FALSE)

MICPath <- paste0(datadir, '/', MICFile)
MICdata1 <- read.table(MICPath, sep = '\t', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
MICdata2 <- MICdata1[,1]
names(MICdata2) <- rownames(MICdata1)
MIC <- MICdata2[order(MICdata2)]

mnames1 <- names(MIC)
mnames <- gsub("_foldInd$", "", mnames1)
names(MIC) <- mnames

cnames1 <- names(expdata1)
cnames <- gsub("_foldInd$", "", cnames1)
names(expdata1) <- cnames



expdata2 <- expdata1[, names(MIC)]
expdata3 <- lapply(expdata2, function(x) as.numeric(x))
cdata1 <- do.call('cbind', expdata3)
rownames(cdata1) <- rownames(expdata2)

rnames1 <- rownames(cdata1)
rnames <- gsub(".*?(A1S_[[:digit:]]+)$", "\\1", rnames1)
cdata <- cdata1
rownames(cdata) <- rnames



lclass1 <- rep("", length(MIC))
names(lclass1) <- names(MIC)
lclass1[MIC <= miccoundary] = "S"
lclass1[MIC > miccoundary] = "R"
lclass <- factor(lclass1)
species = 'Acinetobacter baumannii'

getAnno <- function(feature) {
    fval <- entrez_search(db="protein", feature)
    for (lid in fval$ids) {
        fval2 <- entrez_summary(db="protein", id = lid)
        ltitle <- fval2$title
        if (grepl(species, ltitle)) {
            res <- gsub("(.*?)[[:space:]]*\\[.*$", "\\1", fval2$title)
            return (paste0(feature, ":", res))
        }
    }
    return (paste0(feature, ":"))
}

fVals <- lapply(rnames, getAnno)
fMap <- unlist(fVals)
names(fMap) <- rnames 


cdata <- cdata[, !(colnames(cdata) %in% "RB197")]
MIC <- MIC[!(names(MIC) %in% "RB197")]

lclass <- lclass[!(names(lclass) %in% "RB197")]

alldata <- list(lclass = lclass, MIC = MIC, MICMid = MICMid, cdata = cdata, species = species, fMap = fMap)

outpath <- paste0(datadir, '/', outfile)
save(alldata, file = outpath)

