# 1.takes list of driver mutations and subsets vcf file
# 2.also counts non-synonymous variants and outputs to table
# 3.inputs summary list of drivers and assigns driver mutation catagory

################### notes ###################


################# subroutines #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterLynchList.WGS.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/sampleList.Set.09.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[1]])
#sampleNames <- sampleNames[-15]

#load driver gene list
annoDriverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table.05-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)

Filedir <- "1.platypusCalls/somaticTotal/"
FilePrep <- ".snv.annoVar.variant_function.txt"

geneNameCol <- 3


for(currSet in 1:length(sampleNames)){
  paste("##### plotting VAF distribution for ", sampleNames[currSet]," ######", sep="")
  
  purFileFlag <- 1
  
  currentSample <- subset(sampleList, sampleList[1]==sampleNames[currSet])
  nonNorSamples <- currentSample[[2]]
  norName <- nonNorSamples[(currentSample[1,7]+1)]
  nonNorSamples <- nonNorSamples[-(currentSample[1,7]+1)]
  noSamples <- nrow(currentSample)
  
  #get SNVs and timings
  fileName <- paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], FilePrep, sep="")
  dataFile <- read.table(file=fileName, header=FALSE, stringsAsFactors = FALSE, sep="\t")
  
  #get allele frequencies
  for(currSam in 1:noSamples){
    dataFile[ncol(dataFile)+1] <- dataFile[[7+noSamples+currSam]] / dataFile[[7+currSam]]
  }
  
  #dataFile[1] <- NULL
  names(dataFile) <- c("type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(currentSample[[2]], ".NV", sep=""), paste(currentSample[[2]], ".NR", sep=""), currentSample[[2]])
  
  #calculate correction factor from clonal variants
  dataFile[ncol(dataFile)+1] <- 0
  names(dataFile)[ncol(dataFile)] <- "CT"
  dataFile[ncol(dataFile)+1] <- 0
  names(dataFile)[ncol(dataFile)] <- "phyloLoc"
  for(currRow in 1:nrow(dataFile)){
    currBases <- dataFile[currRow, c("ref", "alt")]
    if((as.character(currBases[1])=="C" & as.character(currBases[2])=="T") | (as.character(currBases[1])=="G" & as.character(currBases[2])=="A")){
      dataFile[currRow, "CT"] <- 1
    }
    assessRow <- table(dataFile[currRow, nonNorSamples] > 0)
    if("TRUE" %in% names(assessRow)){
      if(as.numeric(assessRow["TRUE"]) == length(nonNorSamples)){
        #trunk
        dataFile[currRow, "phyloLoc"] <- "T"
      }else if(as.numeric(assessRow["TRUE"]) == 1){
        #leaf
        dataFile[currRow, "phyloLoc"] <- "L"
      }else{
        #branch
        dataFile[currRow, "phyloLoc"] <- "B"
      }
    }
  }  
  dataFile <- dataFile[dataFile[["CT"]]==1, ]
  
  
  #plot VAF distribution of C>T variants
  pdf(file=paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], ".CT.VAF.pdf", sep=""), onefile=TRUE, width=10, height=(noSamples*4))
  par(mfrow=c(noSamples, 2), xpd=TRUE, mar=c(4,4,4,4))
  
  for(j in 1:noSamples){
    tempHolder <- dataFile[, currentSample[j,2]]
    tempHolder <- tempHolder[tempHolder>0]

    hist(tempHolder, breaks=seq(0,1,0.01), main=currentSample[j,2], xlab="VAF", ylab="f")
  }
  dev.off()
  
  
  #plot depth-VAF distribution of C>T variants
  pdf(file=paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], ".CT.depth.pdf", sep=""), onefile=TRUE, width=10, height=(noSamples*4))
  par(mfrow=c(noSamples, 2), xpd=TRUE, mar=c(4,4,4,4))
  
  for(j in 1:noSamples){
    tempHolder <- dataFile[, c(paste(currentSample[j,2], ".NV", sep=""), currentSample[j,2])]
    
    tempHolder <- tempHolder[tempHolder[[2]]>0, ]
    
    plot(x=tempHolder[[1]], y=tempHolder[[2]], main=currentSample[j,2], xlab="depth", ylab="VAF")
  }
  dev.off()
  
  
}  


  
  