# script compared original conformance Platypus calls (non-synonymous) with Platypus bulk recalling
# outputs a comparison summary table with proportions of recalled variants
# additionally reconstructs a conformance mutect table for comparison

#################### notes #################### 


################## libraries ##################


################# subroutines ################# 


################# main program ################# 

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
mutectSampleList <- read.csv(file="~/PhD/CRCproject/masterBulkList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)

SNVdir <- "1.platypusCalls/somaticTotal.0.05/"
mutDir <- "1.platypusCalls/bulkCalls/"
mutOutDir <- "1c.bulkRecall/"

#platApp <- ".snv.annoVar.variant_function.txt"
platApp <- ".snv.annoVar.variant_function.0.05.txt"

#get set names
setNames <- unique(mutectSampleList[[1]])

#make recal table
recallSummary <- data.frame(matrix(NA, ncol=length(setNames), nrow=9))
names(recallSummary) <- setNames
row.names(recallSummary) <- c("totalPlatypus", "presentInBulk", "recalled", "TrunkVar", "BranchVar", "LeafVar", "reTrunkVar", "reBranchVar", "reLeafVar")

#load driver gene list
driverList <- read.csv(file="~/PhD/CRCproject/9.driverEvents/driverList.chromLoci.csv", header=FALSE, stringsAsFactors=FALSE)

for(currSet in 1:length(setNames)){
  
  #subset sampleList
  currSampleList <- subset(sampleList, sampleList[1]==setNames[currSet])
  normIndex <- currSampleList[1,7]+1
  currNames <- currSampleList[[2]]
  noSamples <- currSampleList[1,8]
  
  #get set platypus file
  platName <- paste(currSampleList[1,6], SNVdir, currSampleList[1,1], "/", currSampleList[1,1], platApp, sep="")
  platypusData <- read.table(file=platName, sep="\t", stringsAsFactors=FALSE, header=FALSE)
  platypusData <- platypusData[-1]
  names(platypusData) <- c("region", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(currNames, ".NV", sep=""), paste(currNames, ".NR", sep=""))
  platypusData["pos2"] <- paste(platypusData[["chrom"]], ":", platypusData[["pos"]], sep="")
  
  #get allele frequencies
  for(currVar in 1:noSamples){
    NVcol <- paste(currNames[currVar], ".NV", sep="")
    NRcol <- paste(currNames[currVar], ".NR", sep="")
    platypusData[ncol(platypusData) + 1] <- platypusData[[NRcol]] / platypusData[[NVcol]] 
    names(platypusData)[ncol(platypusData)] <- currNames[currVar]
  }
  
  #get driver gene status, parse gene names to new column
  for(l in 1:nrow(platypusData)){
    geneName <- strsplit(platypusData[l, "gene"], "," )
    platypusData[l,2] <- geneName[[1]][1]
  }
  
  #separate gene names
  for(n in 1:nrow(platypusData)){
    platypusData[n, "gene"] <- strsplit(platypusData[n, "gene"], ":")[[1]][1]
  }
  
  #get driver gene color scheme for plot
  platypusData[ncol(platypusData)+1] <- NA
  names(platypusData)[ncol(platypusData)] <- "driver"
  for(n in 1:nrow(platypusData)){
    if(platypusData[n, "gene"] %in% driverList[[5]]){
      platypusData[n, ncol(platypusData)] <- "driver"
    }
  }
  
  #get non-normal sample list
  nonNorSam <- currNames[-normIndex]
  noSamples <- noSamples-1
  
  #get mutect bulk sample data in list
  mutName <- paste(currSampleList[1,6], mutDir, setNames[currSet], "/", setNames[currSet], ".snv.annoVar.variant_function.txt", sep="")
  mutectData <- read.table(file=mutName, sep="\t", stringsAsFactors=FALSE, header=FALSE)
  mutectData[ncol(mutectData)+1] <- paste(mutectData[["contig"]], mutectData[["position"]], sep=":")
  names(mutectData) <- c("region", "gene", "chrom", "pos", "pos2", "ref", "alt", "N.NR", "B.NR", "N.NV", "B.NV", "lociPos")
  mutectData["lociPos"] <- paste(mutectData[["chrom"]], ":", mutectData[["pos"]], sep="") 
  
  #make recal table
  recallTable <- data.frame(matrix(NA, ncol=2, nrow=nrow(platypusData)))
  names(recallTable) <- c("recalled", "phyloLoc")
  row.names(recallTable) <- platypusData[["pos2"]]

  #for each variant search in mutect files and populate recal table 
  for(currVar in 1:nrow(recallTable)){
    #get platypus data
    platVar <- paste(platypusData[currVar, "chrom"], ":", platypusData[currVar, "pos"], sep="")
    assessRow <- table(platypusData[currVar, nonNorSam] > 0)
    if(as.numeric(assessRow["TRUE"]) == noSamples){
      recallTable[currVar, 2] <- "T"
    }else if(as.numeric(assessRow["TRUE"]) == 1){
      recallTable[currVar, 2] <- "L"
    }else{
      recallTable[currVar, 2] <- "B"
    }
    
    #look up sample
    if(platVar %in% mutectData[["lociPos"]]){
      #variant was found, check allele freq
      if(mutectData[mutectData[["lociPos"]]==platVar, "B.NV"] > 0){
        recallTable[currVar, 1] <- 1
      }else{
        print(paste("variant", platVar,"in set", setNames[currSet], "has zero allele frequency"))
      }
    }else{
      #variant was present in sample but not found in bulk
      recallTable[currVar, 1] <- 0
    }
  }
    
  #store recall data
  recallSummary[1, currSet] <- nrow(platypusData)
  recallSummary[2, currSet] <- nrow(recallTable[recallTable[["recalled"]]==1, ])
  recallSummary[3, currSet] <- recallSummary[2, currSet] / recallSummary[1, currSet]
  recallSummary[4, currSet] <- table(recallTable[[2]])["T"]
  recallSummary[5, currSet] <- table(recallTable[[2]])["B"]
  recallSummary[6, currSet] <- table(recallTable[[2]])["L"]
  recallSummary[7, currSet] <- nrow(recallTable[recallTable[["phyloLoc"]]=="T" & recallTable[["recalled"]]==1, ]) / recallSummary[4, currSet]
  recallSummary[8, currSet] <- nrow(recallTable[recallTable[["phyloLoc"]]=="B" & recallTable[["recalled"]]==1, ]) / recallSummary[5, currSet]
  recallSummary[9, currSet] <- nrow(recallTable[recallTable[["phyloLoc"]]=="L" & recallTable[["recalled"]]==1, ]) / recallSummary[6, currSet]
}
  
#write final table  
outSumName <- paste(currSampleList[1,6], mutOutDir, "totalRecallSummary.platypus.csv", sep="")
write.table(recallSummary, file=outSumName, sep=",", row.names=TRUE, col.names=TRUE, quote=FALSE)  


