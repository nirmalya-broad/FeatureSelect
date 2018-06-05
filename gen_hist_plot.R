library(ggpubr)
library(ggplot2)
tab_file <- '/home/nirmalya/research/DataDx/kpmero_random.csv'
pdf_file <- '/home/nirmalya/research/DataDx/kpmero_random.pdf'

tab_file <- '/home/nirmalya/research/DataDx/kpcip_random.csv'
pdf_file <- '/home/nirmalya/research/DataDx/kpcip_random.pdf'
ltab <- read.table(tab_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
temp = data.frame(table(cut(x = ltab$accuracy, breaks = seq(0.5, 1, 0.025))))
p <- ggplot(ltab) + geom_histogram(aes(x = accuracy), breaks = seq(0.5,1,0.025), 
    fill = replace(rep("grey", NROW(temp)), which.max(temp$Freq), "black")) +
    labs(x = "Accuracy", y = "Frequency") + theme(axis.text = element_text(size = 13), axis.title = element_text(size = 13))

ggsave(pdf_file, p)


