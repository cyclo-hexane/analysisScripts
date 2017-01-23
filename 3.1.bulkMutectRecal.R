# script compared original conformance Platypus calls (non-synonymous) with Mutect bulk recalling
# outputs a comparison summary table with proportions of recalled variants
# additionally reconstructs a conformance mutect table for comparison

#################### notes #################### 


################## libraries ##################


################# subroutines ################# 


################# main program ################# 

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
mutectSampleList <- read.csv(file="~/PhD/CRCproject/1c.bulkRecall/masterSampleList.bulks.csv", header=FALSE, stringsAsFactors=FALSE)

SNVdir <- "1.platypusCalls/nonSyn.0.05/"
mutDir <- "1b.mutectCalls/"
mutOutDir <- "1c.bulkRecall/"

#platApp <- ".snv.annoVar.variant_function.txt"
platApp <- ".snv.annoVar.exonic_variant_function.0.05.txt"

#get set names
setNames <- mutectSampleList[[1]]

#make recal table
recallSummary <- data.frame(matrix(NA, ncol=length(setNames), nrow=12))
names(recallSummary) <- setNames
row.names(recallSummary) <- c("totalPlatypus", "presentInBulk", "recalled", "TrunkVar", "BranchVar", "LeafVar", "reTrunkVar", "reBranchVar", "reLeafVar", "trunkDiff", "branchDiff", "LeafDiff")

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
  names(platypusData) <- c("region", "gene", "chrom", "pos", "pos2", "ref", "alt", currNames, paste(currNames, ".NV", sep=""))
  platypusData["pos2"] <- paste(platypusData[["chrom"]], ":", platypusData[["pos"]], sep="")
  
  #get allele frequencies
  for(currVar in 1:noSamples){
    NVcol <- paste(currNames[currVar], ".NV", sep="")
    platypusData[ncol(platypusData) + 1] <- platypusData[[NVcol]] / platypusData[[currNames[currVar]]]
    names(platypusData)[ncol(platypusData)] <- paste(currNames[currVar], ".AF", sep="")
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
  mutName <- paste(currSampleList[1,6], mutDir, setNames[currSet], "/", mutectSampleList[currSet, 3], sep="")
  mutectData <- read.table(file=mutName, sep="\t", stringsAsFactors=FALSE, header=TRUE)
  mutectData[ncol(mutectData)+1] <- paste(mutectData[["contig"]], mutectData[["position"]], sep=":")
  names(mutectData)[ncol(mutectData)] <- "lociPos"
   
  
  #make recal table
  recallTable <- data.frame(matrix(NA, ncol=noSamples, nrow=nrow(platypusData)))
  names(recallTable) <- nonNorSam
  row.names(recallTable) <- platypusData[["pos2"]]
  platAFcol <- paste(nonNorSam, ".AF", sep="")
  
  #for each variant search in mutect files and populate recal table 
  for(currVar in 1:nrow(recallTable)){
    #get platypus data
    platVar <- platypusData[currVar, platAFcol]
    
    #look in each mutect file for each sample
    for(currLook in 1:noSamples){
      if(platVar[1, currLook] == 0){
        #mark on table as absent in platypus data
        recallTable[currVar, currLook] <- 0
        next
      }else{
        #look up sample
        if(row.names(recallTable[currVar,]) %in% mutectData[["lociPos"]]){
          #variant was found
          recallTable[currVar, currLook] <- 2
        }else{
          #variant was present in sample but not found in bulk
          recallTable[currVar, currLook] <- 1
        }
      }
    }
  }
    
  #store recall data
  totalVariants <- table(unlist(recallTable))
  recallSummary[1, currSet] <- totalVariants["1"] + totalVariants["2"]
  recallSummary[2, currSet] <- totalVariants["2"]
  recallSummary[3, currSet] <- totalVariants["2"] / recallSummary[1, currSet]
  
  #assess phylogenetic change
  recallTable[ncol(recallTable) + 1] <- NA
  names(recallTable)[ncol(recallTable)] <- "phyloLoc"
  for(currPhylo in 1:nrow(recallTable)){
    assessRow <- table(recallTable[currPhylo, ] > 0)
    if(assessRow["TRUE"] == 1){
      recallTable[currPhylo, "phyloLoc"] <- "L"
    }else if(assessRow["TRUE"] == noSamples){
      recallTable[currPhylo, "phyloLoc"] <- "T"
    }else{
      recallTable[currPhylo, "phyloLoc"] <- "B"
    }
  }
  phyloTab <- table(recallTable["phyloLoc"])
  recallSummary[4, currSet] <- phyloTab["T"]
  recallSummary[5, currSet] <- phyloTab["B"]
  recallSummary[6, currSet] <- phyloTab["L"]
  
  #output table
  outName <- paste(currSampleList[1,6], mutOutDir, currSampleList[1,1], "/",currSampleList[1,1], ".bulkRecal.txt", sep="")
  write.table(recallTable, file=outName, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
  
  #output table
  outNameDrive <- paste(currSampleList[1,6], mutOutDir, currSampleList[1,1], "/",currSampleList[1,1], ".driverRecal.txt", sep="")
  driverTable <- recallTable
  driverTable[ncol(driverTable)+1] <- platypusData["driver"]
  driverTable[ncol(driverTable)+1] <- platypusData["gene"]
  driverTable <- driverTable[!is.na(driverTable[["driver"]]), ]
  write.table(driverTable, file=outNameDrive, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
   
  #now assess for mutect recaptured variants
  recallTable[ncol(recallTable)] <- NA
  for(currPhylo in 1:nrow(recallTable)){
    assessRow <- table(recallTable[currPhylo, ] == 2)
    if("TRUE" %in% names(assessRow)){
      if(assessRow["TRUE"] == noSamples){
        recallTable[currPhylo, "phyloLoc"] <- "T"
      }else if(assessRow["TRUE"] == 1){
        recallTable[currPhylo, "phyloLoc"] <- "L"
      }else{
        recallTable[currPhylo, "phyloLoc"] <- "B"
      }
    }else{
      recallTable[currPhylo, "phyloLoc"] <- "notFound"
    }
  }
  
  phyloTab <- table(recallTable["phyloLoc"])
  recallSummary[7, currSet] <- phyloTab["T"]
  recallSummary[8, currSet] <- phyloTab["B"]
  recallSummary[9, currSet] <- phyloTab["L"]
  
  #get proportional change
  recallSummary[10, currSet] <- phyloTab["T"] / recallSummary[4, currSet]
  recallSummary[11, currSet] <- phyloTab["B"] / recallSummary[5, currSet]
  recallSummary[12, currSet] <- phyloTab["L"] / recallSummary[6, currSet]
  
}
  
#write final table  
outSumName <- paste(currSampleList[1,6], mutOutDir, "/totalRecallSummary.csv", sep="")
write.table(recallSummary, file=outSumName, sep=",", row.names=TRUE, col.names=TRUE, quote=FALSE)  


