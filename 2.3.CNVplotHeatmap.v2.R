# uses CNV disruptions from cloneHD segmentations

##########################   notes   ##########################
#
#   1.populate CNV loci and gene list  
#   
#   2.order table by lesions
#   
#   3.plot summary of CNVs         
#

########################## subroutines  ##########################

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

#### get p and q arm locations from list
#chrPos <- subChromPos; centPos <- centromerePos
getPQ <- function(chrPos, centPos){
  pEnd <- max(chrPos[chrPos[[3]] <= centPos[1], "end"])
  qStart <- min(chrPos[chrPos[[3]] > centPos[2], "start"])
  posVect <- c(0, pEnd, qStart, chrPos[nrow(chrPos), "end"])
  return(posVect)
}

#get copy state for a given CNV
#stateData <- subSegCNVs; sampleSubNames <- noNorNames   
getCNVstate <- function(stateData, sampleSubNames){
  subStateList <- c()
  for(currTest in 1:length(sampleSubNames)){
    tempColNames <- c(paste(sampleSubNames[currTest], "_Minor", sep=""), paste(sampleSubNames[currTest], "_Major", sep=""))
    if((stateData[1, tempColNames[1]] == 1 & stateData[1, tempColNames[2]] == 2) | (stateData[1, tempColNames[1]] == 0 & (stateData[1, tempColNames[1]] + stateData[1, tempColNames[2]]) == 3)){
      subStateList[currTest] <- 3
    }else if((stateData[1, tempColNames[1]] == 1 & stateData[1, tempColNames[2]] > 1) | (stateData[1, tempColNames[1]] + stateData[1, tempColNames[2]]) >= 4){
      subStateList[currTest] <- 4
    }else if(stateData[1, tempColNames[1]] == 0 & stateData[1, tempColNames[2]] == 1){
      subStateList[currTest] <- 1
    }else if(stateData[1, tempColNames[1]] == 0 & stateData[1, tempColNames[2]] == 2){
      subStateList[currTest] <- 2
    }else{
      subStateList[currTest] <- 0
    }
  }
  
  #determine the most common CNV state
  subStateList <- subStateList[subStateList != 0]
  subStateList <- table(subStateList)
  maxStates <- as.numeric(names(subStateList)[as.numeric(subStateList) == max(as.numeric(subStateList))])
  
  #if states are equally prepresented take the highest value
  if(length(maxStates) > 1){
    maxStates <- max(maxStates)
  }
  
  return(maxStates)
}

######################### libraries ##########################

######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[[1]])
setNames <- setNames[c(1:20)]

#working dir
CNVholdingDir <- "7.CNVcalls.final/baf.updated/"
CNVfileNames <- ".penalty0.95.baf.gt.txt"

#chromosome positions
chromPos <- read.table(file="~/PhD/ReferenceGenome/cytoBand.txt", stringsAsFactors = FALSE, header=FALSE, sep="\t")
names(chromPos) <- c("chrom", "start", "end", "desc")
chromPos[4] <- paste(chromPos[[1]], chromPos[[4]], sep="")

#get table to annotate
cnvTab <- read.csv(file="~/PhD/CRCproject/9.driverEvents/CNVs/CNVfreqTable.csv", stringsAsFactors = FALSE, header = TRUE)
#cnvTab <- cnvTab[cnvTab[[2]]=="N", ]

plotDir <- "9.driverEvents/CNVs/"


##################### 1.populate CNV table ##################### 
for(currSeg in 1:length(setNames)){
  print(paste("#### getting segs for sample ", setNames[currSeg], " ####",sep=""))
  setName <- setNames[currSeg]

  #set names for all samples in current set
  currentSamples <- sampleList[sampleList[[1]]==setName, ]
  samNames <- currentSamples[[2]]
  noSamples <- length(samNames)-1
  noNorNames <- samNames[-(currentSamples[1,7]+1)]
  
  if(currentSamples[1,1] == "Ha.4" | currentSamples[1,1] == "Ha.6" | currentSamples[1,1] == "Ha.12" | currentSamples[1,1] == "Ha.22"){
    noNorNames <- paste("X", noNorNames, sep = "")
  }
  
  #setup input/output names
  segCNVfile <- paste(sampleList[currSeg, 6], CNVholdingDir, setName, CNVfileNames, sep="")
  segCNV <- read.table(file=segCNVfile, sep="\t", fill=TRUE, stringsAsFactors=FALSE, header=TRUE)
  
  #set ambigous calls to WT
  segCNV[segCNV=="X"] <- 1
  
  #save sampled genome length then filter small regions
  segCNV["nloci"] <- (segCNV[["last.locus"]] - segCNV[["first.locus"]]) / 1000000
  segCNV <- segCNV[segCNV[["nloci"]] > 1, ]
  
  #ensure states are numeric
  for(currChg in 5:ncol(segCNV)){
    segCNV[currChg] <- as.numeric(segCNV[[currChg]])
  }
  
  ##### A. filter CNVs for valid regions and assign phylogenetic location
  returnedInfo <- filterCNVs(segCNV, noSamples)
  filteredSegCNVs <- returnedInfo[[1]]
  filteredSegCNVs <- filteredSegCNVs[filteredSegCNVs[["phyloLoc"]]=="T", ]
  
  #assess disruptions across each entry on cnvTab
  for(currEnt in 1:nrow(cnvTab)){

    #get chromosome name and arm coordinates
    splitChrom <- cnvTab[currEnt, 1]
    if(length(which(strsplit(splitChrom, "")[[1]]=="p")) != 0){
      splitChrom <- strsplit(splitChrom, split = "p")[[1]][1]
      subChromPos <- chromPos[chromPos[[1]]==splitChrom, ]
      centromerePos <- subChromPos[subChromPos[[5]]=="acen", ]
      subChromPos <- subChromPos[subChromPos[["end"]] <= centromerePos[1, "start"], ]
    }else{
      splitChrom <- strsplit(splitChrom, split = "q")[[1]][1]
      subChromPos <- chromPos[chromPos[[1]]==splitChrom, ]
      centromerePos <- subChromPos[subChromPos[[5]]=="acen", ]
      subChromPos <- subChromPos[subChromPos[["start"]] > centromerePos[1, "end"], ]
    }
    currChr <- splitChrom
    testCoord <- c(subChromPos[1, "start"], subChromPos[nrow(subChromPos), "end"])

    currLociSize <- (testCoord[2] - testCoord[1]) / 1000000
    
    #subset CNVs for current chromosome arm
    subSegCNVs <- filteredSegCNVs[filteredSegCNVs[["chr"]]==currChr, ]
    if(nrow(subSegCNVs) == 0){
      next
    }
    subSegCNVs <- subSegCNVs[subSegCNVs[["first.locus"]] >= testCoord[1] & subSegCNVs[["last.locus"]] <= testCoord[2], ]
    if(nrow(subSegCNVs) == 0){
      next
    }
    
    
    #determine % of current arm disrupted
    sumDisruption <- sum(subSegCNVs[["nloci"]])
    
    #mark as arm disruption if greater than 50% of the arm is effected
    if((sumDisruption / currLociSize) > 0.4){
      #get the majority CNV state
      currState <- getCNVstate(subSegCNVs, noNorNames)
      cnvTab[currEnt, (currSeg+2)] <- currState
    }
  
  }
}



##################### 2.order events on table ##################
cnvTab[ncol(cnvTab)+1] <- NA
cnvTab[ncol(cnvTab)+1] <- NA
cnvTab[ncol(cnvTab)+1] <- NA
cnvTab[ncol(cnvTab)+1] <- NA
names(cnvTab)[(ncol(cnvTab)-3):ncol(cnvTab)] <- c("loss", "LOH", "Tri", "Tet")
for(currApp in 1:nrow(cnvTab)){
  currAss <- table(as.numeric(cnvTab[currApp, setNames]))
  currAss <- currAss[as.numeric(names(currAss)) > 0]
  if(1 %in% as.numeric(names(currAss))){
    cnvTab[currApp, "loss"] <- as.numeric(currAss["1"])
  }
  if(2 %in% as.numeric(names(currAss))){
    cnvTab[currApp, "loss"] <- as.numeric(currAss["2"])
  }
  if(3 %in% as.numeric(names(currAss))){
    cnvTab[currApp, "Tri"] <- as.numeric(currAss["3"])
  }
  if(4 %in% as.numeric(names(currAss))){
    cnvTab[currApp, "Tet"] <- as.numeric(currAss["4"])
  }
}
cnvTab[ncol(cnvTab)+1] <- NA
for(x in 1:nrow(cnvTab)){
  if(sum(cnvTab[x, c("Tri", "Tet")], na.rm = TRUE) > 0 & sum(cnvTab[x, c("LOH", "loss")], na.rm = TRUE) == 0){
    cnvTab[x, ncol(cnvTab)] <- 1
  }else if(sum(cnvTab[x, c("Tri", "Tet")], na.rm = TRUE) == 0 & sum(cnvTab[x, c("LOH", "loss")], na.rm = TRUE) > 0){
    cnvTab[x, ncol(cnvTab)] <- 2
  }else{
    cnvTab[x, ncol(cnvTab)] <- 0
  }
}
cnvTab <- cnvTab[order(cnvTab[[ncol(cnvTab)]]), ]
cnvTab <- cnvTab[, -c((ncol(cnvTab)-4):ncol(cnvTab))]

#order by diploidy status
colOrder <- c(4,6,7,2,5,10,1,8,3,9,11,12,15,14,19,17,20,16,18,13)
colOrder <- colOrder + 2
cnvTab <- cnvTab[, c(1,2,colOrder)]

##################### 2.plot summary of CNVs #####################
outputLoc <- paste(sampleList[1,6], plotDir, "CNVlociFreq.pdf", sep="")
pdf(file=outputLoc, width=5, height=6)
  par(mar=c(2,5,8,2), xpd=TRUE, cex.main = 2.5) 
  plot(1, 1, col="white", axes=F, xlim=c(0, ncol(cnvTab)), ylim=c(0,nrow(cnvTab)), xlab="", ylab="", main="")
  
  #plot segmentations
  for(currSam in 1:(ncol(cnvTab)-2)){
    for(chromo in 1:nrow(cnvTab)){
      #assign box colour
      if(cnvTab[chromo, (currSam+2)] == 0){
        plotCol <- "grey"
      }else if(cnvTab[chromo, (currSam+2)] == 2){
        plotCol <- "blue"
      }else if(cnvTab[chromo, (currSam+2)] == 1){
        plotCol <- "steelblue"
      }else if(cnvTab[chromo, (currSam+2)] == 3){
        plotCol <- "red"
      }else{
        plotCol <- "red4"
      }
      
      #plot square 
      rect(xleft=currSam, xright=(currSam+1), ytop=chromo, ybottom=(chromo+1), col=plotCol)
      
    }
  }
  
  axis(side=3, at=seq(1,(ncol(cnvTab)-2), 1), labels=names(cnvTab)[3:ncol(cnvTab)], las=2, lty = 0)
  axis(side=2, at=c(1:nrow(cnvTab)), labels=cnvTab[[1]], las=2, lty = 0, line=-2)
dev.off()


