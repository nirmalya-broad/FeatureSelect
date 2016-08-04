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
infile <- 'KpMero1_norm_foldInd.txt'
MICFile <- 'KpMero1_strains_MICs.txt'
outfile <- 'KpMero.RData'
MICMid <- 1

inpath <- paste0(datadir, '/', infile)
expdata1 <- read.csv(inpath, header = TRUE, sep = " ", row.names = 1, stringsAsFactors = FALSE)

MICPath <- paste0(datadir, '/', MICFile)
MICdata1 <- read.table(MICPath, sep = ' ', header = TRUE, row.names = 1, stringsAsFactors = FALSE)
MICdata2 <- MICdata1[,2]
names(MICdata2) <- MICdata1[,1]
MIC <- MICdata2[order(MICdata2)]
lclass_d <- MICdata1[,3]
names(lclass_d) <- MICdata1[,1]
lclass1 <- lclass_d
lclass1[lclass1 == "I"] <- "R"

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
rnames <- gsub(".*?(KPN_[[:digit:]]+).*$", "\\1", rnames1)
cdata <- cdata1
rownames(cdata) <- rnames



#lclass1 <- rep("", length(MIC))
#names(lclass1) <- names(MIC)
#lclass1[MIC <= miccoundary] = "S"
#lclass1[MIC > miccoundary] = "R"
lclass2 <- factor(lclass1)

lclass <- lclass2[names(MIC)]

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
save(alldata, file = outpath)

