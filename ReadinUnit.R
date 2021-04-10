library(dplyr)

List = c('MG211', 'MG214', 'MG222', 'MG223', 'MG337', 'MG479', 'MG521', 'MG526', 'MG561', 'MG735', 'MG755', 'MG852', 'MG1050', 'MG1053', 'MG1215', 'MG2129', 'MG2158', 'MG2210', 'MG2214', 'MG2321', 'MG2358', 'MG2371', 'MG2382', 'MG2861', 'MG3278', 'MG3295', 'MG3309', 'MG3337', 'MG3346', 'MG3347', 'MG3390', 'MG3994', 'MG4025', 'MG4031', 'MG4063', 'MG4268', 'MG4287', 'MG4494', 'MG4553', 'MG4570', 'MG4932', 'MG5153', 'MG5233', 'MG5316', 'MG5550')


NameList <- readLines(paste( "./inputFiles/IGC_Full_Estimate1/211/save/MG94_01_02_nonclockqq_summary.txt", sep = ""))
NameList <- unlist(strsplit(NameList[90], split=" "))
NameList <- NameList[!(NameList %in% NameList[1])]

NameList4Swap <- readLines(swapSamplePath)
NameList4Swap <- unlist(strsplit(NameList4Swap[count], split = " "))
NameList4Swap <- NameList4Swap[!(NameList4Swap %in% NameList4Swap[1])]



Readin <- function(list, iterations, caseName) {
  
  columnNames = c()
  data = c()
  for (i in List) {
    
    folderName = paste("/", substr(i, start = 3, stop = nchar(i)), sep = "")
    
    path = paste("./inputFiles/", caseName, folderName, "/save/MG94_01_02_nonclockqq_summary.txt", sep = "")
    numbers <- NA
    tryCatch(numbers <- readLines(path,  n = 89) %>% as.numeric(), error=function(e){
      message(" No data yet.")
    })
    
    variableName = paste(i, iterations, sep = "")
    columnNames <- append(columnNames, variableName)
    data <- cbind(data, numbers)
    
  }
  data <- matrix(data = data, ncol = 45, nrow = 89, byrow = FALSE)
  data.frame(data)
  colnames(data) <- columnNames
  rownames(data) <- NameList
  return(data)
}



ReadinF <- function(list, iterations, caseName) {
  
  columnNames = c()
  data = c()
  for (i in List) {
    
    folderName = paste("/", substr(i, start = 3, stop = nchar(i)), sep = "")
    
    path = paste("./inputFiles/", caseName, folderName, "/save/Force_MG94_01_02_nonclockqq_summary.txt", sep = "")
    numbers <- NA
    tryCatch(numbers <- readLines(path,  n = 89) %>% as.numeric(), error=function(e){
      message(" No data yet.")
    })
    
    variableName = paste(i, iterations, sep = "")
    columnNames <- append(columnNames, variableName)
    data <- cbind(data, numbers)
    
  }
  data <- matrix(data = data, ncol = 45, nrow = 89, byrow = FALSE)
  data.frame(data)
  colnames(data) <- columnNames
  rownames(data) <- NameList
  return(data)
}


count = 0
swapSamplePath = 'inputFiles/NewUnswappedClocked3/211/save/MG94_01_02_clock_summary.txt'
for (i in readLines(swapSamplePath)) {
  count = count + 1
}




Readin4Swap <- function(list, iterations, caseName) {
  
  columnNames = c()
  data = c()
  for (i in List) {
    
    folderName = paste("/", substr(i, start = 3, stop = nchar(i)), sep = "")
    
    path = paste("./inputFiles/", caseName, folderName, "/save/MG94_01_02_nonclockqq_summary.txt", sep = "")
    numbers <- NA
    tryCatch(numbers <- readLines(path,  n = (count - 1)) %>% as.numeric(), error=function(e){
      message(" No data yet.")
    })
    
    variableName = paste(i, iterations, sep = "")
    columnNames <- append(columnNames, variableName)
    data <- cbind(data, numbers)
    
  }
  data <- matrix(data = data, ncol = 45, nrow = (count - 1), byrow = FALSE)
  data.frame(data)
  colnames(data) <- columnNames
  rownames(data) <- NameList4Swap
  return(data)
}




Readin4SwapClocked <- function(list, iterations, caseName) {
  
  columnNames = c()
  data = c()
  for (i in List) {
    
    folderName = paste("/", substr(i, start = 3, stop = nchar(i)), sep = "")
    
    path = paste("./inputFiles/", caseName, folderName, "/save/MG94_01_02_clock_summary.txt", sep = "")
    numbers <- NA
    tryCatch(numbers <- readLines(path,  n = (count - 1)) %>% as.numeric(), error=function(e){
      message(" No data yet.")
    })
    
    variableName = paste(i, iterations, sep = "")
    columnNames <- append(columnNames, variableName)
    data <- cbind(data, numbers)
    
  }
  data <- matrix(data = data, ncol = 45, nrow = (count - 1), byrow = FALSE)
  data.frame(data)
  colnames(data) <- columnNames
  rownames(data) <- NameList4Swap
  return(data)
}




Readin4SwapFClocked <- function(list, iterations, caseName) {
  
  columnNames = c()
  data = c()
  for (i in List) {
    
    folderName = paste("/", substr(i, start = 3, stop = nchar(i)), sep = "")
    
    path = paste("./inputFiles/", caseName, folderName, "/save/Force_MG94_01_02_clock_summary.txt", sep = "")
    numbers <- NA
    tryCatch(numbers <- readLines(path,  n = (count - 1)) %>% as.numeric(), error=function(e){
      message(" No data yet.")
    })
    
    variableName = paste(i, iterations, sep = "")
    columnNames <- append(columnNames, variableName)
    data <- cbind(data, numbers)
    
  }
  data <- matrix(data = data, ncol = 45, nrow = (count - 1), byrow = FALSE)
  data.frame(data)
  colnames(data) <- columnNames
  rownames(data) <- NameList4Swap
  return(data)
}

