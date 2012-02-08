
readPWMfile <- function(file)
{
  table <- read.table(file=file, fill=T)
  table.vector <- cbind(as.vector(table[[1]]), as.vector(table[[2]]), as.vector(table[[3]]), as.vector(table[[4]]), as.vector(table[[5]]))

  XX <- which(table.vector=="XX")
  DE <- which(table.vector=="DE")
  names <- table.vector[DE,2]
  listPWM <- list()

  for (i in seq(length(DE)))
  {
    listPWM[[i]] <- matrix(as.numeric(table.vector[(DE[i]+1):(XX[i]-1),2:5]), nrow=4, byrow=T, dimnames=list(c("A","C","G","T")))
    colnames(listPWM[[i]]) <- 1:(length(listPWM[[i]])/4)
    names(listPWM)[i] <- table.vector[DE[i],2]
  }
  return (listPWM)
}


