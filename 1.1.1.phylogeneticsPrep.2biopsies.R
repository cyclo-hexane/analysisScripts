# 1. takes a list of 2 biopsy samples and outputs lengths of branches in a 'Y' tree

################### notes ###################


################# subroutines #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterAdinCaList.exc.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[21:24]

#load driver gene list
annoDriverList <- read.csv(file="~/PhD/cryptProject/driverGeneRef.csv", header=TRUE, stringsAsFactors=FALSE)

Filedir <- "1.platypusCalls/somaticTotal/"
outDir <- "2.phylogenetics/somaticTotal.0.05/"
FilePrep <- ".snv.annoVar.variant_function.txt"
outPrep <- ".phyloLocs.txt"

for(currSet in 1:length(sampleNames)){
  currentSample <- subset(sampleList, sampleList[1]==sampleNames[currSet])
  nonNorSamples <- currentSample[[2]]
  nonNorSamples <- nonNorSamples[-(currentSample[1,7]+1)]
  noSamples <- length(nonNorSamples)
  
  #get variants
  fileName <- paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], FilePrep, sep="")
  dataFile <- read.table(file=fileName, header=FALSE, stringsAsFactors = FALSE, sep="\t")
  names(dataFile) <- c("type", "gene", "chrom", "pos", "pos1", "ref", "alt", paste(currentSample[[2]], ".NR", sep=""), currentSample[[2]])
  
  branchLengths <- matrix(0, nrow=1, ncol=3)
  colnames(branchLengths) <- c("trunk", nonNorSamples[1], nonNorSamples[2])
  
  #loop through variants and count up branch length
  for(currRow in 1:nrow(dataFile)){
    if(dataFile[currRow, nonNorSamples[1]] > 0 & dataFile[currRow, nonNorSamples[2]] > 0){
      branchLengths[1, "trunk"] <- branchLengths[1, "trunk"] + 1
    }else if(dataFile[currRow, nonNorSamples[1]] > 0){
      branchLengths[1, nonNorSamples[1]] <- branchLengths[1, nonNorSamples[1]] + 1
    }else{
      branchLengths[1, nonNorSamples[2]] <- branchLengths[1, nonNorSamples[2]] + 1
    }
  }
  
  #write total table
  phyloFileOut <- paste(currentSample[1,6], outDir, currentSample[1,1], "/", currentSample[1,1], outPrep, sep="")
  write.table(branchLengths, file=phyloFileOut, sep="\t", quote = FALSE, row.names=FALSE, col.names = TRUE)
  
}



