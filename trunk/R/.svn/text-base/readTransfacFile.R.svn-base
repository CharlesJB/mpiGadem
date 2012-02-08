
readTransfacFile <- function(x)
{
	jasp <- read.table(file=x, fill=T)
	jaspar <- cbind(as.vector(jasp[[1]]), as.vector(jasp[[2]]), as.vector(jasp[[3]]), as.vector(jasp[[4]]), as.vector(jasp[[5]]))

	XX <- which(jaspar=="XX")
	DE <- which(jaspar=="DE")
	names <- jaspar[DE,2]
	PWM <- list()

	for (i in 1:length(DE))
	{
		PWM[[i]] <- matrix(as.numeric(jaspar[(DE[i]+1):(XX[i]-1),2:5]), nrow=4, byrow=T, dimnames=list(c("A","C","G","T")))
		colnames(PWM[[i]]) <- 1:(length(PWM[[i]])/4)
		names(PWM)[i] <- jaspar[DE[i],2]
	}
return (PWM)
}






