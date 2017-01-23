# 1.takes list of driver mutations and subsets vcf file
# 2.also counts non-synonymous variants and outputs to table
# 3.inputs summary list of drivers and assigns driver mutation catagory

################### notes ###################


################# subroutines #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/sampleList.Set.09.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[1]])
#sampleNames <- sampleNames[-15]

#load driver gene list
annoDriverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table.05-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)

Filedir <- "1.platypusCalls/somaticTotal/"
#Filedir <- "1.platypusCalls/exome/"
FilePrep <- ".snv.annoVar.variant_function.txt"
#FilePrep <- ".snv.annoVar.exonic_variant_function.txt"
outPrep <- ".snv.Platypus.clonal.txt"

geneNameCol <- 2

SNVcountTable <- data.frame(matrix(NA, nrow=nrow(sampleList), ncol=2))
SNVcountTable[1] <- sampleList[[2]]
names(SNVcountTable) <- c("biopsyName", "total.variants")

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
  #dataFile <- dataFile[-1]
  
  #get allele frequencies
  for(currSam in 1:noSamples){
    dataFile[ncol(dataFile)+1] <- dataFile[[7+noSamples+currSam]] / dataFile[[7+currSam]]
  }
  
  #dataFile[1] <- NULL
  names(dataFile) <- c("type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(currentSample[[2]], ".NV", sep=""), paste(currentSample[[2]], ".NR", sep=""), currentSample[[2]])
  
  dataFile <- dataFile[c("type", "gene", "chrom", "pos", "pos2", "ref", "alt", currentSample[[2]] )]
  
  #count variants and add to list
  for(currBio in nonNorSamples){
    SNVcountTable[SNVcountTable[[1]]==currBio, 2] <- nrow(dataFile[dataFile[[currBio]] > 0,])
  }
 
  #calculate correction factor from clonal variants
  dataFile[ncol(dataFile)+1] <- NA
  names(dataFile)[ncol(dataFile)] <- "phyloLoc"
  for(currRow in 1:nrow(dataFile)){
   assessRow <- table(dataFile[currRow, currentSample[[2]]] > 0)
   if("TRUE" %in% names(assessRow)){
     if(assessRow["TRUE"] == length(nonNorSamples)){
       #trunk
       dataFile[currRow, "phyloLoc"] <- "T"
     }else if(assessRow["TRUE"] == 1){
       #leaf
       dataFile[currRow, "phyloLoc"] <- "L"
     }else{
       #branch
       dataFile[currRow, "phyloLoc"] <- "B"
     }
    }
  }
  #subset for truncal variants
  trunkVar <- dataFile[dataFile[["phyloLoc"]]=="T", ]
  
  #get purity estimates
  purName <- paste(currentSample[1,6], "13.purityEstimates/", currentSample[1,1], ".penalty0.95.purity-cna-baf.txt", sep="")
  if(file.exists(purName)){
    purData <- read.table(file=purName, header=FALSE, stringsAsFactors = FALSE, sep=" ")
    
    #overwrite with means based estiamtes
    for(currOver in 1:nrow(purData)){
      currT <- purData[currOver, 3]
      purData[currOver, 4] <- mean(trunkVar[[currT]])*2
    }
  }else{
    purFileFlag <- 0
  }
  
  
  pdf(file=paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], ".VAF.pdf", sep=""), onefile=TRUE, width=10, height=(noSamples*4))
  par(mfrow=c(noSamples, 2), xpd=TRUE, mar=c(4,4,4,4))
  
  for(j in 1:noSamples){
    tempHolder <- dataFile[, currentSample[j,2]]
    tempHolder <- tempHolder[tempHolder>0]
    
    #current purity correction
    if(currentSample[j,2] == norName | purFileFlag == 0){
      purCorr <- 0
    }else{
      purCorr <- purData[purData[[3]]==currentSample[j,2], 4]
      purCorr <- 1/purCorr
    }
    tempHolderCorr <- purCorr * tempHolder
    tempHolderCorr[tempHolderCorr > 1] <- 1

    hist(tempHolder, breaks=seq(0,1,0.01), main=currentSample[j,2], xlab="VAF", ylab="f")
    hist(tempHolderCorr, breaks=seq(0,1,0.01), main=paste("corrected ", currentSample[j,2], sep=""), xlab="VAF", ylab="f")
  }
  dev.off()
  
  
  #plot just truncal variants
  pdf(file=paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], ".truncal.VAF.pdf", sep=""), onefile=TRUE, width=10, height=(noSamples*4))
  par(mfrow=c(noSamples, 2), xpd=TRUE, mar=c(4,4,4,4))
  
  for(j in 1:noSamples){
    tempHolder <- trunkVar[, currentSample[j,2]]
    tempHolder <- tempHolder[tempHolder>0]
    
    #current purity correction
    if(currentSample[j,2] == norName | purFileFlag == 0){
      purCorr <- 0
    }else{
      purCorr <- purData[purData[[3]]==currentSample[j,2], 4]
      purCorr <- 1/purCorr
    }
    tempHolderCorr <- purCorr * tempHolder
    tempHolderCorr[tempHolderCorr > 1] <- 1
    
    hist(tempHolder, breaks=seq(0,1,0.01), main=currentSample[j,2], xlab="VAF", ylab="f")
    hist(tempHolderCorr, breaks=seq(0,1,0.01), main=paste("corrected ", currentSample[j,2], sep=""), xlab="VAF", ylab="f")
  }
  dev.off()
  
}  
 

#plot depth against VAF
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
    dataFile[ncol(dataFile)+1] <- dataFile[[8+noSamples+currSam]] / dataFile[[8+currSam]]
  }
  
  dataFile[1] <- NULL
  names(dataFile) <- c("type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(currentSample[[2]], ".NV", sep=""), paste(currentSample[[2]], ".NR", sep=""), currentSample[[2]])
  
  #dataFile <- dataFile[c("type", "gene", "chrom", "pos", "pos2", "ref", "alt", currentSample[[2]] )]
  
  #calculate correction factor from clonal variants
  dataFile[ncol(dataFile)+1] <- NA
  names(dataFile)[ncol(dataFile)] <- "phyloLoc"
  for(currRow in 1:nrow(dataFile)){
    assessRow <- table(dataFile[currRow, currentSample[[2]]] > 0)
    if("TRUE" %in% names(assessRow)){
      if(assessRow["TRUE"] == length(nonNorSamples)){
        #trunk
        dataFile[currRow, "phyloLoc"] <- "T"
      }else if(assessRow["TRUE"] == 1){
        #leaf
        dataFile[currRow, "phyloLoc"] <- "L"
      }else{
        #branch
        dataFile[currRow, "phyloLoc"] <- "B"
      }
    }
  }
  #subset for truncal variants
  trunkVar <- dataFile[dataFile[["phyloLoc"]]=="T", ]

  
  pdf(file=paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], ".depth-VAF.pdf", sep=""), onefile=TRUE, width=10, height=(noSamples*4))
  par(mfrow=c(noSamples, 2), xpd=TRUE, mar=c(4,4,4,4))
  
  for(j in 1:noSamples){
    tempHolder <- dataFile[, c(paste(currentSample[j,2], ".NV", sep=""), currentSample[j,2])]
    tempHolder <- tempHolder[tempHolder[[2]]>0, ]
    
    plot(tempHolder[[1]], tempHolder[[2]], main=currentSample[j,2], xlab="depth", ylab="VAF")
  }
  dev.off()
  
  
  #plot just truncal variants
  pdf(file=paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], ".truncal.depth-VAF.pdf", sep=""), onefile=TRUE, width=10, height=(noSamples*4))
  par(mfrow=c(noSamples, 2), xpd=TRUE, mar=c(4,4,4,4))
  
  for(j in 1:noSamples){
    tempHolder <- trunkVar[, c(paste(currentSample[j,2], ".NV", sep=""), currentSample[j,2])]
    tempHolder <- tempHolder[tempHolder[[2]]>0, ]
    
    plot(tempHolder[[1]], tempHolder[[2]], main=currentSample[j,2], xlab="depth", ylab="VAF")
  }
  dev.off()
  
}  

SNVcountName <- paste(currentSample[1,6], Filedir, "SNVcounts.txt", sep="")
#SNVcountName <- paste(currentSample[1,6], Filedir, "SNVcounts.exome.txt", sep="")
write.table(SNVcountTable, file=SNVcountName, sep="\t", row.names = FALSE, quote = FALSE)

  