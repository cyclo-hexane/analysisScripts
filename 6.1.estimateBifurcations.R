# estimates trunk bifurcation from variant call and CNV sets

################### notes ###################
#
#   6.0.CNVtvaraintClustering.R (est alpha/beta and make input table: X.matlab.input.txt)
#             |
#             V
#   AnalyseCNVData.m (matlab script on apocrita, output: X.matlab.output.txt)
#             |
#             V
#   6.1.estimateBifurcations.R (used for plot)
#             |
#             V
#   6.2.assessCNVtiming.heatmap.R (make final figure of timings)
#
################# libraries #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
setNames <- c("Set.02", "Set.03", "Set.04", "Set.05", "Set.09.Distal", "Set.09.Proximal", "Set.10", "Polyp.08.WGS")

CNVphyloDir <- "7.CNVcalls/baf.updated/"
phyloPrep <- ".penalty0.95.baf.gt.txt"

SNVDir <- "1.platypusCalls/somaticTotal.0.05/"
SNVPrep <- ".snv.annoVar.variant_function.0.05.txt"

#all timings table
bifurcationEst <- data.frame(matrix(NA, ncol=length(setNames), nrow=3))
names(bifurcationEst) <- setNames
rownames(bifurcationEst) <- c("length", "totalVariants", "truncalVariants")

#table output
biDir <- "6.CNVtiming/inputData/"

#loop through sample sets
for(currSam in 1:length(setNames)){
  print(paste("##### getting bifurcation for ", setNames[currSam], " #####", sep=""))
  
  #get SNV calls
  currentNames <- subset(sampleList, sampleList[1]==setNames[currSam])
  sampleNames <- currentNames[[2]]
  
  noSamples <- currentNames[1,8]
  normalIndex <- currentNames[1,7]+1
  nonNormalNames <- currentNames[-normalIndex, 2]
  
  #get SNV file
  SNVfileName <- paste(currentNames[1,6], SNVDir, setNames[currSam], "/", setNames[currSam], SNVPrep, sep="")
  confData <- read.table(file=SNVfileName, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  names(confData) <- c("no", "region", "gene", "chrom", "pos", "pos2", "alt", "ref", paste(sampleNames, ".NR", sep=""), paste(sampleNames, ".NV", sep=""))
  
  #get CNV calls
  CNVfileName <- paste(currentNames[1,6], CNVphyloDir, setNames[currSam], phyloPrep, sep="")
  CNVdata <- read.table(file=CNVfileName, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
  #get diploid regions within CNV calls
  keepRegions <- c()
  keepCounter <- 1
  for(currCNV in 1:nrow(CNVdata)){
    assessCNV <- CNVdata[currCNV, c(5:ncol(CNVdata))]
    assTab <- table(assessCNV==1)
    if("TRUE" %in% names(assTab)){
      if(as.numeric(assTab["TRUE"]) == ((noSamples-1)*2)){
        keepRegions[keepCounter] <- currCNV
        keepCounter <- keepCounter + 1
      }
    }
  }
  CNVdata <- CNVdata[keepRegions, ]
  CNVdata["nloci"] <- CNVdata[["last.locus"]] - CNVdata[["first.locus"]]
  
  #if no regions are diploid, estimation cannot be made
  if(nrow(CNVdata)==0){
    next
  }
  
  #calculate allele frequencies
  for(y in 1:noSamples){
    confData[(ncol(confData)+1)] <- as.vector(confData[(8+noSamples+y)]) / as.vector(confData[8+y])
    names(confData)[ncol(confData)] <- paste(currentNames[y, 2], ".AF", sep="")
  }
  
  #assess phylo location
  confData[ncol(confData)+1] <- NA
  names(confData)[ncol(confData)] <- "phyloLoc"
  for(currRow in 1:nrow(confData)){
    currAssess <- confData[currRow, paste(nonNormalNames, ".AF", sep="")]
    currNo <- table(currAssess!=0)
    if("TRUE" %in% names(currNo)){
      if(as.numeric(currNo["TRUE"]) == (noSamples-1)){
        confData[currRow, "phyloLoc"] <- "T"
      }else{
        confData[currRow, "phyloLoc"] <- "N"
      }
    }
  }

  #remove non-diploid regions from truncal variants
  keepSNVs <- c()
  keepCounter <- 1
  for(currRow in 1:nrow(confData)){
    currChrom <- confData[currRow, "chrom"]
    currPos <- confData[currRow, "pos"]
    
    #subset CNV diploid list
    subCNVS <- subset(CNVdata, CNVdata[1]==currChrom)
    
    if(nrow(subCNVS)>=1){
      #check to see if SNV in CNV regions
      for(currRegion in 1:nrow(subCNVS)){
        currStart <- subCNVS[currRegion, "first.locus"]
        currEnd <- subCNVS[currRegion, "last.locus"]
        
        if(currPos > currStart & currPos < currEnd){
          keepSNVs[keepCounter] <- currRow
          keepCounter <- keepCounter + 1
          next
        }
      }
    }
  }
  confData <- confData[keepSNVs, ]
  
  #now store values
  bifurcationEst[1, currSam] <- sum(as.numeric(CNVdata[["nloci"]]))
  bifurcationEst[2, currSam] <- nrow(confData)
  bifurcationEst[3, currSam] <- nrow(confData[confData[["phyloLoc"]]=="T", ])
}

#save file
biFurcationFile <- paste(currentNames[1,6], biDir, "bifurcationTable.txt", sep="")
write.table(bifurcationEst, file=biFurcationFile, sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
