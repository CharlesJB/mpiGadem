setMethod("[",
"gadem",
function(x, i, j=ANY, ..., exact=TRUE, drop=FALSE){
    selected <- NULL
    selectedMatch <- NULL
    if (is.numeric(i))
    {
    selected = i
    } else {
        if (exact)
        {
            selected = which(names(x) %in% i)
        } else {
            selected = unlist(sapply(i, function(x){grep(x, names(gadem))}))
        }    
    }
    if (drop)
    {selected <- (1:length(x))[-selected]}
    for (l in unique(selected))
    {
        selectedMatch <- c(selectedMatch, new("motif", alignList=x@motifList[[l]]@alignList ,consensus=x@motifList[[l]]@consensus, pwm=x@motifList[[l]]@pwm, name=x@motifList[[l]]@name))          
    }        
    if (!is.null(selectedMatch))
    {        res <- new("gadem", motifList=selectedMatch,parameters=x@parameters)
    } else {
        res <- NULL
    }                
    return(res)    
})

setMethod("[[","gadem",
    function(x, i, j, ..., exact = TRUE)

    {
      if(length(i) != 1)
      {
        stop("subscript out of bounds (index must have length 1)")
      }
      if(missing(i))
      {
        return(x)
      }
      if(!missing(j))
      {
        stop("incorrect number of dimensions")
      }
      x@motifList[[i]]
})


setMethod("dim",
"gadem",
function(x){
  .Deprecated("nOccurrences", package="rGADEM")
sapply(x@motifList,function(x){length(x@alignList)})
})

setGeneric("nOccurrences",function(x) standardGeneric("nOccurrences"))
setMethod("nOccurrences",
"gadem",
function(x){
sapply(x@motifList,function(x){length(x@alignList)})
})

setGeneric("nMotifs",function(x) standardGeneric("nMotifs"))
setMethod("nMotifs",
"gadem",
function(x){
length(x@motifList)
})

setMethod("names",
"gadem",
function(x){
  sapply(x@motifList,function(x){x@name})
})


setGeneric("consensus",function(x) standardGeneric("consensus"))
setMethod("consensus",
"gadem",
function(x){
  sapply(x@motifList,function(x){x@consensus})
})

setGeneric("parameters",function(x) standardGeneric("parameters"))
setMethod("parameters",
"gadem",
function(x){
x@parameters
})


setGeneric("getPWM",function(x) standardGeneric("getPWM"))
setMethod("getPWM",
"motif",
function(x){
  pwm<-x@pwm
return(pwm)
})

setMethod("getPWM",
"gadem",
function(x){
  pwm<-lapply(x@motifList,"getPWM")
  names(pwm)<-names(x)
return(pwm)
})

setMethod("plot",
"motif",
function(x,y="MISSING",...){
pwm<-makePWM(getPWM(x))
plot(pwm,...)
})

setMethod("plot",
"gadem",
function(x,y="MISSING",...){
  x<-lapply(x@motifList,function(x,...){plot(makePWM(x@pwm),...)},...)
})


setGeneric("startPos", function(x) standardGeneric("startPos"))

setMethod("startPos",
"gadem",
function(x){
    
 	ListStart <- list()
    for (i in seq(length(x)))
    {
	ListStart[[i]]=x@motifList[[i]]@alignList[[i]]@start
	names(ListStart)[i]=names(x)[i]      
     }        
   
    return(ListStart)
})


setGeneric("endPos", function(x) standardGeneric("endPos"))
setMethod("endPos",
"gadem",
function(x){
    
 	ListEnd <- list()
    for (i in seq(length(x)))
    {
	ListEnd[[i]]=x@motifList[[i]]@alignList[[i]]@end
	names(ListEnd)[i]=names(x)[i]      
     }        
   
    return(ListEnd)
})

###############################
############SHOW###############
###############################

#####motiv#####
setMethod("show", "gadem",
function(object)
{
	cat("\tObject of class 'gadem'","\n")
	cat("\tThis object has the following slots: \n")
	cat("\tmotifs,pwm,consensus,align,name,seq,chr,start,end,strand,seqID,pos,pval,fastaHeader\n\n")
})










