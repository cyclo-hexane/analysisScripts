# script filters SNV file for CNV regions


################### notes ###################

#libraries


################### subroutines ####################



################### main program ####################

#sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)

#holdingDir <- arguments[2]
holdingDir <- "2.phylogenetics/somaticTotal.0.05/"

#vcfName <- arguments[3]
varName <- ".snv.annoVar.variant_function.0.05.txt"

#CNV file
CNVholding <- "8.CNVphylogenetics/CNVphylo/"

CNVname <- ".phyloCNVs.csv"

outDir <- "2.phylogenetics/somaticTotalCNVfiltered/"
outName <- ".snv.annoVar.variant_function.0.05.txt"

setNames <- unique(sampleList[[1]])

for(currSam in 1:length(setNames)){
  subsample <- subset(sampleList, sampleList[1]==setNames[currSam])
  
  #get SNV file
  SNVfileName <- paste(subsample[1,6], holdingDir, subsample[1,1], "/", subsample[1,1], varName, sep="")
  SNVin <- read.table(file=SNVfileName, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  
  #get CNV file
  CNVfileName <- paste(subsample[1,6], CNVholding, subsample[1,1], "/", subsample[1,1], CNVname, sep="")
  CNVin <- read.csv(file=CNVfileName, stringsAsFactors = FALSE, header=TRUE)
  
  deleteVect <- c()
  deleteCounter <- 1
  for(currRow in 1:nrow(SNVin)){
    currCNVs <- subset(CNVin, CNVin["chr"]==SNVin[currRow, 4])
    if(nrow(currCNVs)==0){
      next
    }
    flag <- 0
    for(cnv in 1:nrow(currCNVs)){
      if(SNVin[currRow, 5] >= currCNVs[cnv, "first.locus"] & SNVin[currRow, 5] <= currCNVs[cnv, "last.locus"] & flag==0){
        #remove from SNV list
        deleteVect[deleteCounter] <- currRow
        deleteCounter <- deleteCounter + 1
        flag <- 1
      }
    }
  }
  
  #remove SNVs in CNV region
  SNVin <- SNVin[-deleteVect, ]
  
  #write new file for phylogenetic construction
  outFileName <- paste(subsample[1,6], outDir, subsample[1,1], "/", subsample[1,1], outName, sep="")
  write.table(SNVin, file=outFileName, sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
}