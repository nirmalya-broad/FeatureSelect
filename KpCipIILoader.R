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
infile <- 'KpCip2_allRespNormAbx_forNirmalya.csv'
MICFile <- 'KpCip2_MIC_table_for_Nirmalya.csv'
outfile <- 'KpCipII.RData'
MICMid <- 0.5

inpath <- paste0(datadir, '/', infile)
expdata1 <- read.csv(inpath, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

MICPath <- paste0(datadir, '/', MICFile)

MICdata1 <- read.table(MICPath, sep = '\t', header = TRUE, stringsAsFactors = FALSE)

MICdata2 <- MICdata1[,2]
names(MICdata2) <- MICdata1[,1]
MIC <- MICdata2[order(MICdata2)]
lclass_d <- MICdata1[,3]
names(lclass_d) <- MICdata1[,1]
lclass1 <- lclass_d
lclass1[lclass1 == "I"] <- "R"

cnames1 <- names(expdata1)
cnames <- gsub("_cip$", "", cnames1)
names(expdata1) <- cnames

expdata2 <- expdata1[, names(MIC)]
expdata3 <- lapply(expdata2, function(x) as.numeric(x))
cdata <- do.call('cbind', expdata3)
rownames(cdata) <- rownames(expdata2)

lclass2 <- factor(lclass1)

lclass <- lclass2[names(MIC)]

probe_names <- rownames(cdata)

species = 'Klebsiella pneumoniae'

getAnno <- function(feature_full) {
    feature <- sub('^\\S+?-(\\S+)_\\S+', '\\1', feature_full)
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

dname = "KpCipII"

alldata <- list(lclass = lclass, MIC = MIC, MICMid = MICMid, cdata = cdata, species = species, fMap = fMap, dname = dname)

outpath <- paste0(datadir, '/', outfile)
save(alldata, file = outpath)

