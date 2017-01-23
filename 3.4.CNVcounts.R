# 4.2.ploidySummary.R   
#
#
#
################### notes ###################
#
# summerizes ploidy across biopsies
#
#
################# libraries #################


################# subroutines #################
#### A. filter valid CNV events and assign phylogenetic location ####
#segDataSub <- segCNV;noSamSub <- noSamples
filterCNVs <- function(segDataSub, noSamSub){
  returnList <- as.list(NA)
  
  #assess segmentations if normal, remove
  removeSegs <- c()
  counterSub <- 1
  segDataSub <- as.data.frame(append(segDataSub, x = NA))
  names(segDataSub)[1] <- "phyloLoc"
  for(currSeg in 1:nrow(segDataSub)){
    assessSeg <- table(segDataSub[currSeg, 6:(5+(noSamSub*2))] == 1)
    if("TRUE" %in% names(assessSeg)){
      if(as.integer(assessSeg["TRUE"])==(2*noSamSub)){
        #don't keep, its diploid
        removeSegs[counterSub] <- currSeg
        counterSub <- counterSub + 1
      }else{
        #keep
      }
    }
  }
  #subset to remove normal regions
  if(length(removeSegs)>0){
    segDataSub <- segDataSub[-removeSegs,]
  }
  
  #assess phylogenetic location
  for(currSeg in 1:nrow(segDataSub)){
    #get current CNV seg into table
    assessVector <- segDataSub[currSeg, 6:(5+(noSamSub*2))]
    assessTable <- data.frame(matrix(NA, ncol=2, nrow=noSamSub, byrow = TRUE))
    seqAp <- seq(1, (noSamSub*2), 2)
    for(currAp in 1:nrow(assessTable)){
      assessTable[currAp, 1] <- assessVector[seqAp[currAp]]
      assessTable[currAp, 2] <- assessVector[seqAp[currAp]+1]
    }
    assessTable[3] <- paste(assessTable[[1]], ":", assessTable[[2]], sep="")
    names(assessTable) <- c("minor", "major", "merged")
    
    #assess phylo by uniqueness of states
    stateTest <- table(assessTable[["merged"]])
    noStates <- length(unique(names(stateTest)))
    
    #assess number of states
    if(noStates == 1){
      #trunkal
      segDataSub[currSeg, 1] <- "T"
    }else if("1:1" %in% names(stateTest)){
      if(stateTest["1:1"] == (noSamSub-1)){
        #true leaf event
        segDataSub[currSeg, 1] <- "L"
      }else if(noStates == 2){
        #its a true branch event
        segDataSub[currSeg, 1] <- "B"
      }else{
        #its a bifuricated branch (two branch events in same region)
        segDataSub[currSeg, 1] <- "BiB"
      }
    }else{
      #if no samples have WT state they are ether branch or trunk
      if(noStates == 2){
        #its a bifuricated trunk (one branch events in same region as a trunk)
        segDataSub[currSeg, 1] <- "BiT"
      }else{
        #phylogenetic loc cannot be determined
        segDataSub[currSeg, 1] <- "unknown >2 bifurcations"
      }
    }
  }
  returnList[[1]] <- segDataSub
  returnList[[2]] <- removeSegs
  return(returnList)
}


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

setNames <- unique(sampleList[[1]])

#out directory
CNVholdingDir <- "7.CNVcalls.final/baf.updated/"
CNVfileNames <- ".penalty0.95.baf.gt.txt"

#output table
CNVcounts <- data.frame(matrix(NA, nrow = length(setNames), ncol = 2))
names(CNVcounts) <- c("names", "CNVcounts")

for(currSam in 1:length(setNames)){
  #subset main list
  subSample <- subset(sampleList, sampleList[1]==setNames[currSam])
  
  subSample <- subSample[-(subSample[1,7]+1), ]
  sampleNames <- subSample[[2]]
  noSamples <- length(sampleNames)
  CNVcounts[currSam, 1] <- subSample[1,1]
  
  #get seg file
  fileName <- paste(subSample[1,6], CNVholdingDir, subSample[1,1], CNVfileNames, sep="")
  if(file.exists(fileName)){
    segCNV <- read.table(file=fileName, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  }else{
    CNVcounts[currSam, 2] <- NA
    next
  }
  
  segFiltered <- filterCNVs(segCNV, noSamples)
  segFiltered <- segFiltered[[1]]
  segFiltered["nloci"] <- (segFiltered[["last.locus"]] - segFiltered[["first.locus"]]) / 1000000
  
  #filter for size
  segFiltered <- segFiltered[segFiltered[["nloci"]] > 1, ]
  CNVcounts[currSam, 2] <- nrow(segFiltered)
  
}

#write table
CNVcountName <- "~/PhD/CRCproject/0.writeup/suppTables/CNVcountStats.txt"
write.table(CNVcounts, file=CNVcountName, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

wilcox.test(CNVcounts[1:10, 2], CNVcounts[c(12:16,20:23), 2])

