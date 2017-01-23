# using .tre files with different AF cutoffs plot homplasy vs tree length
# from the HI information tree balance is assessed against homoplasy

#### notes ####
#intervals are currently fixed to seq(0,0.26,0.02)

###################### libraries ###################### 

library(apTreeshape)
library(ape)

#input sampleList from commandline arguments
arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=3){
	stop("\n#### please use syntax > Rscript 6A.phyloStatsPlot.R < sample list file > < holding directory > < prepended.name > ####\n")
}



###################### subroutines ###################### 

#### get homoplasy value from log file ####
getHI <- function(logFileSub){
	#check if log file exists
	if(file.exists(logFileSub) & file.info(logFileSub)$size > 610){
		HIvalue <- system(command=paste("grep \'HI excluding\' ", logFileSub, sep=""), intern=TRUE)
		if(is.na(HIvalue[1])){
			HIstrings <- NA	
			return(HIstrings)
		}else{
			HIstrings <- strsplit(HIvalue, split=" ")
			return(HIstrings[[1]][6])
		}	
	}else{
		HIstrings <- NA	
		return(HIstrings)
	}
}


#### get tree lengths from log file ####
getLength <- function(logFileSub){
	if(file.exists(logFileSub)){
		HIvalue <- system(command=paste("grep \'Tree length\' ", logFileSub, sep=""), intern=TRUE)
		HIstrings <- strsplit(HIvalue, split=" ")
		return(HIstrings[[1]][4])
	}else{
		HIstrings <- NA
		return(HIstrings)	
	}
}


#### get tree starting lengths from log file ####
getLengthStart <- function(logFileSub){
	if(file.exists(logFileSub)){
		HIvalue <- system(command=paste("grep \'Data matrix has\' ", logFileSub, sep=""), intern=TRUE)
		HIstrings <- strsplit(HIvalue, split=" ")
		return(HIstrings[[1]][6])
	}else{
		HIstrings <- NA
		return(HIstrings)	
	}	
}


#### get balance of tree from multi-phylo file ####
getTreeBalance <- function(treeFileSub, normalName){
	treeList <- as.list(0)
	
	if(file.exists(treeFileSub)){
		#get trees from file
		treeList <- read.nexus(file= treeFileSub)
		
		#no of trees
		notrees <- length(treeList)
		
		if(notrees > 10){
			notrees <- 10
		}
		
		#setup results table
		resultsTable <- as.data.frame(matrix(NA, nrow=notrees, ncol=3))
		names(resultsTable) <- c("tree length", "colless test (yule)", "colless test (PDA)")
		row.names(resultsTable) <- paste("tree", c(1:notrees))
		
		#for each tree get stats
		for(j in 1:notrees){
			
			#convert to single tree object
			phyloTree <- treeList[[j]]
			
			#get tree length
			resultsTable[j,1] <- sum(phyloTree$edge.length)
			
			#drop normal sample
			treeTemp <- drop.tip(phyloTree, normalName)
			
			#convert to treeshape object
			phyloShape <- as.treeshape(treeTemp)
			
			#fix polytomy
			if(length(phyloShape) == 0){
				phyloShape <- as.treeshape(treeTemp, model="pda")
			}
			
			#perform stats
			yuleTemp <- colless.test(phyloShape, model="yule", n.mc=1000, alternative="greater")
			pdaTemp <- colless.test(phyloShape, model="pda", n.mc=1000, alternative="greater")
			
			#populate table
			resultsTable[j,2] <- yuleTemp$p.value
			resultsTable[j,3] <- pdaTemp$p.value
		}
		resultsTable <- resultsTable[ order(resultsTable[1]), ]
		return(resultsTable)	
	}else{
		return(0)
	}
}




###################### main program ###################### 

sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
holdingDir <- arguments[2]
namePrepended <- arguments[3]

#sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
#holdingDir <- "2.phylogenetics/stratified/"
#namePrepended <- ".log"

#get sample names
sampleNames <- unique(sampleList[1])
sampleNamesPub <- c("cancer 1", "cancer 2", "cancer 3", "cancer 4 (MSI)", "cancer 5", "cancer 6", "cancer 7", "cancer 8", "cancer 9 distal", "cancer 9 proximal", "cancer 10", "adenoma 1", "adenoma 2", "adenoma 3", "adenoma 4", "adenoma 4 WGS", "adenoma 5" )

#stats table for fo vs f0.05 comparison
statsTable <- data.frame(matrix(0, nrow=nrow(sampleNames), ncol=2))
names(statsTable) <- c("f0", "f0.06")


#loop through sample sets and analyse trees
for(i in 1:nrow(sampleNames)){
	print(paste("performing analysis for ", sampleNames[i,1], sep=""))
	
  #setup plotting table
  plotTable <- data.frame(matrix(0, ncol=length(seq(0,0.35,0.02)), nrow=100))
  names(plotTable) <- seq(0,0.35,0.02)
  row.names(plotTable) <- c(1:100)
  
	#subset sampleList
	currentSamples <- subset(sampleList, sampleList[1]==sampleNames[i,1])
	
  #for each VAF filter
	counter <- 1
	for(cutoff in seq(0,0.35,0.02)){
		
    #for each bootstrap value, populate table
    for(bootStr in 1:100){
      #get log file name
      logName <- paste(currentSamples[1,6], holdingDir, sampleNames[i,1], "/", sampleNames[i,1], "_f", cutoff, "_", bootStr, namePrepended, sep="")
      
      #add HI value for this AF and bootstrap tree to table
      HIindex <- getHI(logName)
      if(is.na(HIindex)){
        plotTable[bootStr, counter] <- NA
      }else{
        plotTable[bootStr, counter] <- as.numeric(HIindex)
      }
    }
    counter <- counter + 1
	}
  
  #get non-bootstrap table
	nonBootTable <- data.frame(matrix(0, ncol=length(seq(0,0.35,0.02)), nrow=1))
	names(nonBootTable) <- seq(0,0.35,0.02)
	counterNB <- 1
  for(cutoff in seq(0,0.35,0.02)){
    logName <- paste(currentSamples[1,6], holdingDir, sampleNames[i,1], "/", sampleNames[i,1], "_f", cutoff, namePrepended, sep="")
    HIindex <- getHI(logName)
    
    #save result to stats table
    if(cutoff == 0 | cutoff == 0.06){
      colName <- paste("f", cutoff, sep="")
      statsTable[i, colName] <- as.numeric(HIindex)
    }
    
    if(is.na(HIindex)){
      nonBootTable[1, counterNB] <- NA
    }else{
      nonBootTable[1, counterNB] <- as.numeric(HIindex)
    }
    counterNB <- counterNB + 1
  }
  #add zero last column
	nonBootTable[1, (ncol(nonBootTable)+1)] <- nonBootTable[1, ncol(nonBootTable)]
  
  #plot results using boxplot
	pdf(file=paste(sampleList[i,6], holdingDir, currentSamples[1,1],".AFfilt.pdf", sep=""), onefile=TRUE, width=20, height=8)
	par(mar=c(8,8,8,8), cex.lab=1.5, cex.axis = 1.5, cex.main=3, cex.lab = 2)
	  plotMed <- boxplot(plotTable, names=names(plotTable), xlab="VAF filter", ylab="HI", main=paste(" VAF cutoff vs homoplasy index (HI) ", sampleNamesPub[i], sep=""))
	  #plotMed <- plotMed$stats
	  #plotMed <- plotMed[3,]
    segments(x0=c(1:18), x1=c(2:19), y0=as.numeric(nonBootTable[1,c(1:18)]), y1=as.numeric(nonBootTable[1,c(2:19)]), lwd=5, col="red")
	dev.off()  
}		

#get comparison stats
statsTable[1] <- as.numeric(statsTable[[1]])
statsTable[2] <- as.numeric(statsTable[[2]])

statsResults <- wilcox.test(statsTable$f0, statsTable$f0.06, paired = TRUE)
statsStrings <- paste("HI of f0 vs f0.05 groups: p-value = ", statsResults$p.value, sep="")
fileStatsOut <- paste(currentSamples[1,6], holdingDir,"hIfiltStats.txt",sep="")

lapply(statsStrings, write, fileStatsOut, append=FALSE)
write.table(statsTable, file=paste(currentSamples[1,6], holdingDir,"hIfiltTable.txt",sep=""), sep="\t", row.names=FALSE, quote = FALSE)
