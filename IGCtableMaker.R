

FileName = "pillar211.pdf"



file = paste(c("/Users/yixuanyang/Desktop/IGC tables/", FileName), collapse = "")
library(gridExtra)
library(datasets)
setwd("/Users/yixuanyang/Desktop")

Value <- readLines(paste( "/Users/yixuanyang/Desktop/IGC_DATA/MG94_01_02_nonclockqq_summary.txt", sep = ""),n = 89)
Value <- as.numeric(Value)
Name <- readLines(paste( "/Users/yixuanyang/Desktop/IGC_DATA/MG94_01_02_nonclockqq_summary.txt", sep = ""))

Name <- unlist(strsplit(Name[90], split=" "))
Name <- Name[!(Name %in% Name[1])]
table <- data.frame(Name, Value)

pdf(file, height = 25)
grid.table(table)
dev.off()


