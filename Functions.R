library(DT)
library(dplyr)
library(ggplot2)
library(stringr)
# getList function will return a list of folder names in a data folder.
getList <- function(directory){
  cwd <- getwd()
  ListOutput = c()
  for(i in c(1:5588)){
    path = paste(cwd,"/inputfiles/",directory,"/",i,"/save/MG94_01_02_nonclockqq_summary.txt", sep = "")
    file = file.exists(path)
    if(file == TRUE){
      ListOutput = append(ListOutput,paste("MG",as.character(i),sep=""))
    }
  }
  return(ListOutput)
}

summaryFile <- function(folder, subfolder){
  savePath <- paste(folder, "/", subfolder, "/save/", sep = "")
  files <- list.files(savePath)
  if(length(files) == 2){
    for(i in files){
      if(grep("summary",i)){
        return(i)
      }
    }
  }
  else{
    return(FALSE)
  }
}

# IGC2Read is a module of List Function. It reads one data outcome as a vector.
IGC2Read <- function(number, directory, file) {
  cwd <- getwd()
  
  path = paste(cwd, "/inputFiles/",directory,"/", number, "/save/", file, sep = "")
  #file = file.exists(path)
  
  v <- scan(path, what = character())
  splitPoint = (length(v) - 1)/2
  value <- as.numeric(v[1:splitPoint])
  names <- v[(splitPoint + 2):length(v)]
  names(value) <- names
  return(value)
  
}

ListFunction <- function(directory,option) {
  folder = paste("inputFiles/", directory, sep = "")
  cases = as.character(sort(as.numeric(list.files(folder))))
  if(option == "List"){
    output = list()
    for(case in cases) {
      file = summaryFile(paste("inputFiles/", directory, sep = ""),case)
      name = paste("MG", case, sep = "")#case
      number = case
      if(file != FALSE){
        data <- IGC2Read(number, directory, file)
        output[[name]]<- data
      }
      else{
        output[[name]]<- NA
        print(paste(case," do not exist", sep = ""))
      }
    }
    return(output)
  }
  if(option == "Matrix"){
    output = c()
    for(case in cases) {
      file = summaryFile(paste("inputFiles/", directory, sep = ""),case)
      if(file != FALSE){
        name = paste("MG", case, sep = "")
        number = case
        data <- IGC2Read(number, directory, file)
        col <- names(data)
        df <- data.frame(col,data)
        colnames(df) <- c("names", name)
        output = merge(df,output,all=TRUE,sort=FALSE)
        output = merge(df,output,all=TRUE,sort=FALSE)
      }
    }
    row <- output[,1]
    output <- output[,c(2:length(output))]
    rownames(output) <- row
    output <- output[length(output):1] %>% as.matrix()
    return(output) 
  }
}


# Compare
comparell <- function(data1,data2){
  n1 <- names(data1)
  n2 <- names(data2)
  n <- unique(c(n1,n2)) #namelist
  d1 <- c()
  d2 <- c()
  eq <- c()
  output <- list()
  for(i in n){
    ll1 <- data1[i][[1]][2]
    ll2 <- data2[i][[1]][2]
    if(is.na(ll1) == T){
      output <- append(output, data2[i])
      print(str_c("Missing ", i, " in data1"))
    }
    else if(is.na(ll2) == T){
      output <- append(output, data1[i])
      print(str_c("Missing ", i, " in data2"))
    }
    else if(ll1 > ll2){
      output <- append(output, data1[i])
      d1 <- append(d1, i)
    }
    else if(ll1 < ll2){
      output <- append(output, data2[i])
      d2 <- append(d2, i)
    }
    else if(ll1 == ll2){
      output <- append(output, data1[i])
      eq <- append(eq, i)
    }
  }
  print("From data1")
  print(d1)
  print("From data2")
  print(d2)
  print("Equal")
  print(eq)
  return(output)
}




