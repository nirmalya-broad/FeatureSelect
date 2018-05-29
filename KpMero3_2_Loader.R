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
infile <- 'KpMero3_runs1-4_allRespNormLog_2.tsv'
outfile <- 'KpMero3_2.RData'
MICMid <- 1

inpath <- paste0(datadir, '/', infile)
expdata1 <- read.csv(inpath, header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

MIC0 <- expdata1[1,]
lclass0 <- expdata1[2,]
strain0 <- expdata1[3,]
expdata2 <- expdata1[-c(1:3),]

lclass1 <- as.character(lclass0)
MIC1 <- as.numeric(MIC0)
strain1 <- as.character(strain0)
strain2 <-gsub("_foldInd$", "", strain1) 
colnames(expdata2) <- strain2

names(lclass1) <- strain2

lclass2 <- factor(lclass1)
names(MIC1) <- strain2
MIC <- sort(MIC1)
lclass <- lclass2[names(MIC)] 

cdata <- data.matrix(expdata2[, names(lclass)])

probe_names <- rownames(cdata)

species = 'Klebsiella pneumoniae'

getAnno <- function(feature_full) {
    feature <- sub('^\\S+?_\\S+?_(\\S+?_\\d+)_\\S+', '\\1', feature_full)
    fval <- entrez_search(db="protein", feature)
    for (lid in fval$ids) {
        fval2 <- entrez_summary(db="protein", id = lid)
        ltitle <- fval2$title
        if (grepl(species, ltitle)) {
            res <- gsub("(.*?)[[:space:]]*\\[.*$", "\\1", fval2$title)
            return (paste0(feature_full, ":", res))
        }
    }
    return (paste0(feature, ":"))
}

fVals <- lapply(probe_names, getAnno)
fMap <- unlist(fVals)
names(fMap) <- probe_names

dname = "KpMero3"

alldata <- list(lclass = lclass, MIC = MIC, MICMid = MICMid, cdata = cdata, species = species, fMap = fMap, dname = dname)

outpath <- paste0(datadir, '/', outfile)
save(alldata, file = outpath)

