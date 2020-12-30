
name = "521original"
FileName = paste(c(name,".pdf"))
file = paste(c("/Users/yixuanyang/Desktop/IGC_DATA", FileName), collapse = "")
library(gridExtra)
library(gtable, grid)
library(datasets)
library(dplyr)
setwd("/Users/yixuanyang/Desktop")

#Value <- readLines(paste( "/Users/yixuanyang/Desktop/IGC_DATA/MG94_01_02_nonclockqq_summary.txt", sep = ""),n = 89)
#Value <- as.numeric(Value)
#Name <- readLines(paste( "/Users/yixuanyang/Desktop/IGC_DATA/MG94_01_02_nonclockqq_summary.txt", sep = ""))

#Name <- unlist(strsplit(Name[90], split=" "))
#Name <- Name[!(Name %in% Name[1])]
#table <- data.frame(Name, Value)
table <- MGc4_i1[,7] %>% as.matrix() %>% round(digits = 4)
table <- data.frame(rownames(table), table)
colnames(table) <- c(name, 'values')
pdf(file, height = 28)
grid.table(table, rows = c(1:89))
dev.off()


#testForSwap  <- readLines("/Users/yixuanyang/Downloads/MG94_01_02_nonclockqq_summary.txt",  n = 89) %>% as.numeric() %>% as.matrix(byrow=TRUE, ncol= 1) 
#rownames(testForSwap) <- NameList

