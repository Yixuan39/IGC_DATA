IGC2Read <- function(number,directory) {
  cwd <- getwd()
  path = paste(cwd, "/inputFiles/",directory,"/", number, "/save/MG94_01_02_nonclockqq_summary.txt", sep = "")
  v <- scan(path, what = character())
  splitPoint = (length(v) - 1)/2
  value <- as.numeric(v[1:splitPoint])
  names <- v[(splitPoint + 2):length(v)]
  names(value) <- names
  return(value)
}

ListFunction <- function(cases, directory) {
  output = list()
  for(case in cases) {
    name = case
    number = substr(case, start=3, stop = nchar(case))
    data <- IGC2Read(number, directory)
    output[[name]] <- data
    
  }
  return(output)
}

IGC2_Full_noforce_noclock <- ListFunction(List2, "IGC2_Full_noforce_noclock")
