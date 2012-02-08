
GADEM<- function (Sequences,seed=1,genome=NULL,verbose=FALSE,numWordGroup=3,numTop3mer=20,numTop4mer=40,numTop5mer=60,numGeneration=5,
	populationSize=100,pValue=0.0002,eValue=0.0,extTrim=1,minSpaceWidth=0,maxSpaceWidth=10,useChIPscore=0,numEM=40,fEM=0.5,widthWt=80,
	fullScan=0,slideWinPWM=6,stopCriterion=1,numBackgSets=10,weightType=0,bFileName="NULL",Spwm=NULL,minSites=-1,maskR=0,nmotifs=25) 
	{

    if(is(Sequences,"RangedData") & is.null(genome))
    {
      stop("You have specified a RangedData object but no genome is specified")
    }

    if(is(Sequences,"RangedData"))
    {
      spSeq<-as.vector(space(Sequences))
      stSeq<-start(Sequences)
      edSeq<-end(Sequences)
      if(verbose)
      {
        cat("Retrieving sequences... ")
      }
      FastaSeq<-getSeq(genome,spSeq,start=stSeq,end=edSeq)
      FastaXstring<-XStringViews(FastaSeq,subjectClass="DNAString")
      if(verbose)
      {
        cat("Done.\n")
      }
    }
    else if(is(Sequences,"XStringViews"))
    {
      FastaXstring<-Sequences
    }
    else if(is(Sequences,"DNAStringSet"))
    {
	  FastaXstring<-XStringViews(Sequences,subjectClass="DNAString")
      FastaXstring<-Sequences	
    }
    else
    {
      stop("Object 'Sequences' should be of type 'XStringViews', 'DNAStringSet' or 'RangedData'")
    }

		FastaSequence<-DNAStringSet(FastaXstring)
		#fastarecords<-XStringSetToFASTArecords(FastaSequence)
		fastaSeqChar<-as.character(FastaSequence)
		fastarecords<-lapply(seq_len(length(fastaSeqChar)), function(i) list(desc=names(fastaSeqChar)[i], seq=fastaSeqChar[[i]]))
		sequenceFasta<-sapply(fastarecords,"tolower")
		accession<-as.integer(1:length(FastaSequence))	

		Lengthfasta<-length(FastaSequence)

    if(verbose)
    {
      cat("*** Start C Programm ***\n")
    }
    
    if(!is.null(seed))
    {
      # Here I save the seed, so that I reset the system at the end
      if(exists(".Random.seed"))
      {
        save.seed <- .Random.seed
      }
      set.seed(seed)
    }
		# Calling C code with .Call
		obj<-.Call("GADEM_Analysis",sequenceFasta,Lengthfasta,accession,as.logical(verbose),numWordGroup,numTop3mer,numTop4mer,numTop5mer,numGeneration,populationSize,
		pValue,eValue,extTrim,minSpaceWidth,maxSpaceWidth,useChIPscore,numEM,fEM,widthWt,fullScan,slideWinPWM,stopCriterion,
		numBackgSets,weightType,bFileName,Spwm,minSites,maskR,nmotifs)

		i<-1
		j<-1
		
		parameter=list()
		parameter[[1]]<-new("parameters",numWordGroup=numWordGroup,numTop3mer=numTop3mer,verbose=as.numeric(verbose),numTop4mer=numTop4mer,numTop5mer=numTop5mer,numGeneration=numGeneration,
		populationSize=populationSize,pValue=pValue,eValue=eValue,extTrim=extTrim,minSpaceWidth=minSpaceWidth,maxSpaceWidth=maxSpaceWidth,useChIPscore=useChIPscore,
		numEM=numEM,fEM=fEM,widthWt=widthWt,fullScan=fullScan,slideWinPWM=slideWinPWM,stopCriterion=stopCriterion,numBackgSets=numBackgSets,weightType=weightType,bFileName=bFileName,
		nSequences=Lengthfasta,maskR=maskR,nmotifs=nmotifs)
		
		list2=list()

		while(i<100&&(!is.null(obj[[i]])))
		{
			list=list()
			length(obj[[i]][[4]][[1]])

			for(j in 1:length(obj[[i]][[4]][[1]]))
			{

				if(is(Sequences,"RangedData"))
				{	
					ind<-as.numeric(obj[[i]][[4]][[5]][[j]])
					list[[j]]<-new("align",seq=obj[[i]][[4]][[1]][[j]],strand=obj[[i]][[4]][[2]][[j]],pos=obj[[i]][[4]][[3]][[j]],pval=obj[[i]][[4]][[4]][[j]],chr=spSeq[ind],start=stSeq[ind],end=edSeq[ind],seqID=obj[[i]][[4]][[6]][[j]],fastaHeader=obj[[i]][[4]][[5]][[j]])
				}

				else if(is(Sequences,"XStringViews"))
				{
					list[[j]]<-new("align",seq=obj[[i]][[4]][[1]][[j]],strand=obj[[i]][[4]][[2]][[j]],pos=obj[[i]][[4]][[3]][[j]],pval=obj[[i]][[4]][[4]][[j]],seqID=obj[[i]][[4]][[6]][[j]],fastaHeader=obj[[i]][[4]][[5]][[j]])			
				}

			}
			matrixPWM<-round(obj[[i]][[2]],4)
			rownames(matrixPWM)<-c("A","C","G","T")
			colnames(matrixPWM)<-seq(dim(matrixPWM)[2])
			list2[[i]]<-new("motif",alignList=list,consensus=obj[[i]][[1]],pwm=matrixPWM,name=obj[[i]][[5]])
			i=i+1
		}
		
		# Reset the seed
		if(exists("save.seed"))
		{
		  .Random.seed<-save.seed
	  }
		
		
		gadem<-new("gadem",motifList=list2,parameters=parameter)
		return(gadem)
	}

