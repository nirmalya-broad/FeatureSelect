library(scales)
library(reshape2)
library(ggplot2)

args = commandArgs(TRUE)

datapath <- '/home/nirmalya/research/DataDx'
codepath <- '/home/nirmalya/research/featureselect'
filename <- args[1]

codeFile <- paste0(codepath, '/FeatureSelect.R')
source(codeFile)

dataFile <- paste0(datapath, '/', filename)
load(dataFile)
dname <- strsplit(filename, "\\.")[[1]][1]

partitionMethod <- 'alternate'
alldata2 <- doPartition(alldata, partitionMethod)

# Add training and test information and also add MIC information to the 
# strain labels.

cdata1 <- alldata2$cdata
cdata <- t(apply(cdata1, 1, rescale))
trainC <- names(alldata2$trainC)
testC <- names(alldata2$testC)
trainC_val <- as.character(alldata2$trainC)
testC_val <- as.character(alldata2$testC)
MIC <- alldata2$MIC
allC <- colnames(cdata)
lclass <- alldata2$lclass
facet_var <- ifelse (lclass == "S", c('Sus'), c('Res'))

trainC_ex <- paste0(trainC, " (FS, ",  MIC[trainC], ")")
names(trainC_ex) <- trainC
testC_ex <- paste0(testC, " (TV, ",  MIC[testC], ")")
names(testC_ex) <- testC

all_ex <- c(trainC_ex, testC_ex)
all_ex2 <- all_ex[allC]

cdata2 <- cdata

cdata3 <- t(cdata2)
cdata4 <- data.frame(cdata3, all_ex2, facet_var)
cdata5 <- cdata4
rownames(cdata5) <- NULL


# We shall now define var_levels in a way; so that the first 10 variables
# would be based on our ranking. The others are based on alphabetical 
# orders. We may even define a facet_grid in between them.

acbmero_list <- c('A1S_2449', 'A1S_0151', 'A1S_1926', 'A1S_1344', 'A1S_1340',
			'A1S_0891', 'A1S_2353', 'A1S_1889', 'A1S_1699', 'A1S_1341')

kpmero_list <- c('KPN_01107', 'KPN_00042', 'KPN_02742', 'KPN_00364', 
				'KPN_00043', 'KPN_04423', 'KPN_03340', 'KPN_04426',
				'KPN_00736', 'KPN_02418')

kpgent_list <- c('KPN_04615', 'KPN_02936', 'KPN_03589', 'KPN_03772',
				'KPN_00663', 'KPN_00411', 'KPN_03031', 'KPN_01000',
				'KPN_03032', 'KPN_01136') 

top_probes <- NULL
if (dname == 'KpMero') {
	top_probes <- kpmero_list
} else if (dname == 'AcbMero') {
	top_probes <- acbmero_list
} else if (dname == 'KpGent') {
	top_probes <- kpgent_list

}	

top_five = top_probes[1:5]
next_five = top_probes[6:10]
var_levels <- sort(rownames(cdata))
probe_diff <- setdiff(var_levels, top_probes)
var_all <- c(top_probes, probe_diff)

cdata6 <- melt(cdata5, id = c("all_ex2", "facet_var"))
cdata6$all_ex2 <- factor(cdata6$all_ex2, levels = all_ex2)
cdata6$variable <- factor(cdata6$variable, levels = var_all)

probe_type = rep("Others", length(cdata6$variable))
probe_type[cdata6$variable %in% top_five] = "Top-five"
probe_type[cdata6$variable %in% next_five] = "Next-five"



cdata7 <- data.frame(cdata6, probe_type)
cdata7$probe_type = factor(cdata7$probe_type, 
		levels = c("Top-five", "Next-five", "Others"))

ltitle <- paste0(dname, ', alternate, rfRFE, Random Forest')

p <- ggplot(cdata7, aes(variable, all_ex2)) + 
	geom_tile(aes(fill = value), colour = "white") +
	scale_fill_gradient(name = "Scaled\nexpression", low = "white", 
		high = "steelblue") + 
	theme(text = element_text(size=14), 
		#axis.text.x = element_text(size=10,angle= 45)) + 
		axis.text.x = element_blank()) + 
	facet_grid(facet_var ~ probe_type, scales = "free", space = "free") +
	labs(x = "Probes", y = "Strains") + 
	ggtitle(ltitle)

lfile <- paste0(dname, '_alternate_rfRFE_rf_heatmap.pdf')
lpath <- paste0(datapath, '/', lfile)
ggsave(lpath, width = 19.3, height = 11)







