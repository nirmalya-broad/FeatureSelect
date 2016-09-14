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
infile <- 'forNirmalya_KpGent1_norm.txt'
#MICFile <- 'forNirmalya_KpGent1_MICs.txt'
MICFile <- 'KpGent1_strain_MIC_table.txt'
outfile <- 'KpGent.RData'
miccoundary <- 2
MICMid <- 2
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

rnames1 <- rownames(expdata2)

rnames <- gsub(".*_(KPN_[[:digit:]]+)_.*", "\\1", rnames1)


rownames(cdata) <- rnames

lclass1 <- rep("", length(MIC))
names(lclass1) <- names(MIC)
lclass1[MIC <= miccoundary] = "S"
lclass1[MIC > miccoundary] = "R"
lclass <- factor(lclass1)

species = 'Klebsiella pneumoniae'
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


alldata <- list(lclass = lclass, MIC = MIC, MICMid = MICMid, cdata = cdata, species = species, fMap = fMap)

outpath <- paste0(datadir, '/', outfile)
save(alldata, file = outfile)
