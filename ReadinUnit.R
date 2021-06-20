library(dplyr)

List = c('MG211', 'MG214', 'MG222', 'MG223', 'MG337', 'MG479', 'MG521', 'MG526', 'MG561', 'MG735', 'MG755', 'MG852', 'MG1050', 'MG1053', 'MG1215', 'MG2129', 'MG2158', 'MG2210', 'MG2214', 'MG2321', 'MG2358', 'MG2371', 'MG2382', 'MG2861', 'MG3278', 'MG3295', 'MG3309', 'MG3337', 'MG3346', 'MG3347', 'MG3390', 'MG3994', 'MG4025', 'MG4031', 'MG4063', 'MG4268', 'MG4287', 'MG4494', 'MG4553', 'MG4570', 'MG4932', 'MG5153', 'MG5233', 'MG5316', 'MG5550')

newNumbers <- c('211', '214', '222', '223', '337', '479', '521', '526', '561', '735', '755', '852', '1050', '1053', '1215', '2129', '2158', '2210', '2214', '2321', '2358', '2371', '2382', '2861', '3278', '3295', '3309', '3337', '3346', '3347', '3390', '3994', '4025', '4031', '4063', '4268', '4287', '4494', '4553', '4570', '4932', '5153', '5233', '5316', '5550') %>% as.matrix() %>% as.numeric()

List2 <- c('MG228', 'MG271', 'MG285', 'MG287', 'MG295', 'MG335', 'MG339', 'MG357', 'MG358', 'MG366', 'MG367', 'MG455', 'MG505', 'MG518', 'MG530', 'MG540', 'MG547', 'MG718', 'MG725', 'MG728', 'MG729', 'MG730', 'MG731', 'MG733', 'MG736', 'MG738', 'MG742', 'MG746', 'MG788', 'MG790', 'MG802', 'MG809', 'MG851', 'MG858', 'MG864', 'MG919', 'MG958', 'MG1056', 'MG1058', 'MG1059', 'MG1060', 'MG1069', 'MG1076', 'MG1079', 'MG1081', 'MG1085', 'MG1086', 'MG1095', 'MG1098', 'MG1115', 'MG1161', 'MG1214', 'MG1238', 'MG1271', 'MG1275', 'MG1276', 'MG1277', 'MG1279', 'MG1280', 'MG1297', 'MG1301', 'MG1309', 'MG1315', 'MG1351', 'MG1363', 'MG1393', 'MG1412', 'MG1431', 'MG1432', 'MG1554', 'MG1557', 'MG2056', 'MG2064', 'MG2065', 'MG2068', 'MG2070', 'MG2088', 'MG2091', 'MG2102', 'MG2116', 'MG2117', 'MG2118', 'MG2127', 'MG2140', 'MG2141', 'MG2142', 'MG2169', 'MG2233', 'MG2337', 'MG2338', 'MG2340', 'MG2379', 'MG2517', 'MG2529', 'MG2535', 'MG2782', 'MG2788', 'MG2794', 'MG2842', 'MG2851', 'MG2852', 'MG2853', 'MG2862', 'MG2863', 'MG2882', 'MG2885', 'MG2904', 'MG2907', 'MG2954', 'MG2957', 'MG2958', 'MG3011', 'MG3014', 'MG3015', 'MG3037', 'MG3286', 'MG3292', 'MG3300', 'MG3308', 'MG3310', 'MG3382', 'MG3388', 'MG3442', 'MG3452', 'MG3453', 'MG3456', 'MG3476', 'MG3486', 'MG3651', 'MG3805', 'MG3814', 'MG3868', 'MG3934', 'MG3950', 'MG3955', 'MG3978', 'MG4002', 'MG4006', 'MG4023', 'MG4050', 'MG4053', 'MG4055', 'MG4060', 'MG4233', 'MG4242', 'MG4290', 'MG4304', 'MG4311', 'MG4312', 'MG4322', 'MG4488', 'MG4499', 'MG4516', 'MG4527', 'MG4568', 'MG4609', 'MG4657', 'MG4686', 'MG4697', 'MG4698', 'MG4700', 'MG4722', 'MG4726', 'MG4883', 'MG4923', 'MG4924', 'MG4927', 'MG4933', 'MG4948', 'MG4978', 'MG4992', 'MG5105', 'MG5112', 'MG5149', 'MG5209', 'MG5232', 'MG5262', 'MG5267', 'MG5308', 'MG5313', 'MG5315', 'MG5317', 'MG5318', 'MG5320', 'MG5325', 'MG5341', 'MG5362', 'MG5370', 'MG5412', 'MG5546', 'MG5552', 'MG5567', 'MG5572', 'MG5583')

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

Readin2 <- function(List2, iterations, caseName) {
  
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

