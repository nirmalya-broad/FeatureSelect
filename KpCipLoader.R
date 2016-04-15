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
infile <- '/KpCip_for_Nirmalya_scaled.txt'
outfile <- 'KpCip.RData'
MICMid <- 2
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

rnames1 <- rownames(cdata1)
rnames <- gsub("(R_up_|R_dn_|C_)(.*)", "\\2", rnames1)
rownames(cdata) <- rnames

species <- 'Klebsiella pneumoniae'

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
