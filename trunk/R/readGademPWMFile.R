
readGademPWMFile <- function (file)
{
	stat1_init <- read.csv(file=file , fill=T, header=F)
	skiplines <- seq(1, dim(stat1_init)[1],5)
	stat1_num <- stat1_init[-skiplines, ]
	stat1_num <- as.character(stat1_num)
	bla <- strsplit(stat1_num,"[\t, ]")

	stat1_names <- as.character(stat1_init[skiplines,])
	stat1_names <- sapply(stat1_names,    function(a){strsplit(a,">")[[1]][2]}    , USE.NAMES=F, simplify=T)

	n=length(stat1_names)

	listPWM <- list()
	for (i in 1:n)
	{
		p=4*(i-1)+1
		listPWM[[i]] <- matrix(as.numeric(c(bla[[p]][-1], bla[[p+1]][-1], bla[[p+2]][-1], bla[[p+3]][-1])) , nrow=4, byrow=T, dimnames=list(c("A","C","G","T")))
		colnames(listPWM[[i]]) <- 1:(length(listPWM[[i]])/4)
		names(listPWM)[i] <- stat1_names[i]
	}
	return(listPWM)
}

