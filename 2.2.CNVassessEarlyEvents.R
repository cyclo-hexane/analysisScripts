# plot summary of CNV disruptions from cloneHD segmentations

##########################   notes   ##########################
#
#   1. reads segmentations summarize in tables
#             |
#             V
#   2. plot clonal vs subclonal gains/losses/LOH
#             
# 
#
########################## libraries  ##########################

######################### subroutines ##########################
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

#function returns cytoband annotation
#subTab <- annoTable[currAnno,]; subBandList <- cytobandList; subDrivers <- annoDriverList
getCytoAnno <- function(subTab, subBandList, subDrivers){
  #subset drivers list
  subsettedDrive <- subDrivers[subDrivers[[1]]==as.numeric(subTab[1, "chr"]), ]
  startEnd <- subTab[1, c("first.locus", "last.locus")]
  
  #subset cyto bands file
  subsetBands <- subBandList[subBandList[[1]]==as.numeric(subTab[1, "chr"]), ]
  subsetBands <- subsetBands[subsetBands[[3]]>as.numeric(startEnd["first.locus"]) & subsetBands[[3]]<=as.numeric(startEnd["last.locus"]), ]
  
  if(nrow(subsetBands)==1){
    #get possible driver genes
    driverGenesSub <- subsettedDrive[subsettedDrive[[2]]>=subsetBands[1,2] & subsettedDrive[[3]]<=subsetBands[1,3], ]
    
    annotation <- subsetBands[1, 4] 
  }else{
    #get possible driver genes
    driverGenesSub <- subsettedDrive[subsettedDrive[[2]]>=subsetBands[1,2] & subsettedDrive[[3]]<=subsetBands[nrow(subsetBands),3], ]
    
    annotation <- paste(subsetBands[1, 4], ":", subsetBands[nrow(subsetBands), 4], sep="") 
  }
  annotationFinal <- as.list(NA)
  annotationFinal[[1]] <- annotation
  if(nrow(driverGenesSub)>0){
    annotationFinal[[2]] <- driverGenesSub
  }else{
    annotationFinal[[2]] <- NA
  }
  
  return(annotationFinal)
}



######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

setNames <- unique(sampleList[[1]])
setNames <- setNames[-c(21:24)]

#working dir
CNVholdingDir <- "7.CNVcalls.final/baf.updated/"
CNVfileNames <- ".penalty0.95.baf.gt.txt"

#out directory
CNVout <- "7.CNVcalls.final/CNVannotation/"

#annotation lists
cytobandList <- read.table(file="~/PhD/ReferenceGenome/cytoBand.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
annoDriverList <- read.csv(file="~/PhD/CRCproject/9.driverEvents/archive/driverList.chromLoci.csv", header=FALSE, stringsAsFactors=FALSE)


#CNV partition cutoff size (Mb), set to number >0 and unhash table to partition data
sizeCutoff <- 5

#plot tables 
plotLargeGains <- data.frame(matrix(NA, ncol=12, nrow=1000))
names(plotLargeGains) <- c("ClonalTet", "BTet", "LTet", "ClonalTri", "BTri", "LTri", "ClonalLOH", "BLOH", "LLOH", "ClonalLoss", "BLoss", "LLoss") 

plotSmallGains <- data.frame(matrix(NA, ncol=12, nrow=1000))
names(plotSmallGains) <- c("ClonalTet", "BTet", "LTet", "ClonalTri", "BTri", "LTri", "ClonalLOH", "BLOH", "LLOH", "ClonalLoss", "BLoss", "LLoss") 


#row locations and counters for plotLargeGains and plotSmallGains tables 
currRowPos <- seq(1, 12, 3)
names(currRowPos) <- c("Tet", "Tri", "LOH", "Loss") 
largeCNVCounter <- c(1,1,1,1,1,1,1,1,1,1,1,1)
names(largeCNVCounter) <- c("ClonalTet", "BTet", "LTet", "ClonalTri", "BTri", "LTri", "ClonalLOH", "BLOH", "LLOH", "ClonalLoss", "BLoss", "LLoss")
smallCNVCounter <- c(1,1,1,1,1,1,1,1,1,1,1,1)
names(smallCNVCounter) <- c("ClonalTet", "BTet", "LTet", "ClonalTri", "BTri", "LTri", "ClonalLOH", "BLOH", "LLOH", "ClonalLoss", "BLoss", "LLoss")


############ 1.make segmentations summary table ###############

for(currSam in 1:length(setNames)){
  #subset main list
  subSample <- subset(sampleList, sampleList[1]==setNames[currSam])
  
  sampleNames <- subSample[[2]]
  sampleNames <- sampleNames[-(subSample[1,7]+1)]
  noSamples <- length(sampleNames)
  
  #get seg file
  segCNV <- read.table(file=paste(subSample[1,6], CNVholdingDir, subSample[1,1], CNVfileNames, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
  #set ambigous calls to WT
  segCNV[segCNV=="X"] <- 1
  
  #save sampled genome length then filter small regions
  segCNV["nloci"] <- (segCNV[["last.locus"]] - segCNV[["first.locus"]]) / 1000000

  #ensure states are numeric
  for(currChg in 5:ncol(segCNV)){
    segCNV[currChg] <- as.numeric(segCNV[[currChg]])
  }
  
  returnedData <- filterCNVs(segCNV, noSamples)
  dataIn <- returnedData[[1]]
  
  #plot segmentation as chromothripsis barchart
  plotChrTab <- data.frame(matrix(NA, nrow=23, ncol=4))
  names(plotChrTab) <- c("breaks", "noStates", "pbreak", "pStates")
  row.names(plotChrTab) <- c(1:23)
  for(currChr in 1:23){
    tempData <- dataIn[dataIn[["chr"]]==currChr, ]
    tempData <- tempData[tempData[["nloci"]] > 5, ]
    
    plotChrTab[currChr, "breaks"] <- nrow(tempData)
    
    #get number of states
    uniqueStates <- c(1)
    for(currSt in 1:length(sampleNames)){
      currState <- paste(sampleNames[currSt], "_Major", sep="")
      currState <- tempData[currState]
      uniqueStates <- unique(c(as.numeric(currState[[1]]), uniqueStates))
    }
    
    plotChrTab[currChr, "noStates"] <- length(uniqueStates)
  }
  #get statistics
  meanStates <- mean(plotChrTab[["noStates"]])
  sdStates <- sd(plotChrTab[["noStates"]])
  meanBreaks <- mean(plotChrTab[["breaks"]])
  sdBreaks <- sd(plotChrTab[["breaks"]])
  for(currBin in 1:nrow(plotChrTab)){
    plotChrTab[currBin, "pStates"] <-  dnorm(x=plotChrTab[currBin, "noStates"], mean = meanStates, sd = sdStates)
    plotChrTab[currBin, "pbreak"] <- dnorm(x=plotChrTab[currBin, "breaks"], mean = meanBreaks, sd = sdBreaks)
  }
  
  
  pdf(file=paste(subSample[1,6], CNVout, "/", subSample[1,1], "chromothripsis.pdf", sep=""), width=5, height=5)
  par(xpd=TRUE, mfrow=c(2,1))
    barplot(plotChrTab[["breaks"]], ylim = c(0,30))
    text(x = seq(2, 28, length.out = 23), y = rep(28, 23), labels = plotChrTab[["pbreak"]], srt=45, cex=0.4)
    
    barplot(-plotChrTab[["noStates"]], ylim=c(-8,0))
    text(x = seq(2, 28, length.out = 23), y = rep(-8, 23), labels = plotChrTab[["pStates"]], srt=45, cex=0.4)
  dev.off()
}
  #copy data to get annotation table
  annoTable <- dataIn
  
  #annotate and rearragen table
  sortOrder <- c("T", "BiT", "B", "BiB", "L", "unknown >2 bifurcations")
  annoTable <- annoTable[annoTable[["nloci"]] >= 5,]
  annoTable <- annoTable[order(match(annoTable[[1]], sortOrder)), ]
  annoTable["cytoBand"] <- NA
  annoTable["driverGenes"] <- NA
  
  #add cytoband
  for(currAnno in 1:nrow(annoTable)){
    annoReturn <- getCytoAnno(annoTable[currAnno,], cytobandList, annoDriverList)
    annoTable[currAnno, "cytoBand"] <- annoReturn[[1]][1]
    
    if(length(annoReturn[[2]]) != 0){
      annoTable[currAnno, "driverGenes"] <- paste(unlist(annoReturn[[2]][5]), collapse = ":")
    }else{
      annoTable[currAnno, "driverGenes"] <- NA
    }
  }
  
  #output annotation table
  #outAnno <- paste(subSample[1,6], CNVout, subSample[1,1], ".CNV.anno.txt", sep="")
  #write.table(annoTable, file=outAnno, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    
  
  #keep only regions with definative trunk, branch or leaf CNVs
  dataIn <- dataIn[dataIn[["phyloLoc"]]=="T" | dataIn[["phyloLoc"]]=="BiT" | dataIn[["phyloLoc"]]=="B" | dataIn[["phyloLoc"]]=="L", ]
  dataIn[(ncol(dataIn)+1):(ncol(dataIn)+5)] <- 0
  names(dataIn)[(ncol(dataIn)-4):ncol(dataIn)] <- c("N", "Tet", "Tri", "LOH", "Loss")
  
  #assess each seg and assign Normal, Tri, Loss, etc counts and phylogenetic location
  for(assessRow in 1:nrow(dataIn)){
    currRow <- dataIn[assessRow, c(paste(sampleNames, "_Major", sep=""), paste(sampleNames, "_Minor", sep=""))]
    
    for(currAss in 1:noSamples){
      currMin <- paste(sampleNames[currAss], "_Minor", sep="")
      currMaj <- paste(sampleNames[currAss], "_Major", sep="")
      if(currRow[1, currMin] == 1 & currRow[1, currMaj] == 1){
        #normal seg in this sample
        dataIn[assessRow, "N"] <- dataIn[assessRow, "N"] + 1
      }else if(currRow[1, currMin] == 0 & currRow[1, currMaj] == 2){
        #BLOH seg in this sample
        dataIn[assessRow, "LOH"] <- dataIn[assessRow, "LOH"] + 1
      }else if(currRow[1, currMin] == 0 & currRow[1, currMaj] == 1){
        #Loss seg in this sample
        dataIn[assessRow, "Loss"] <- dataIn[assessRow, "Loss"] + 1
      }else if(currRow[1, currMin] == 1 & currRow[1, currMaj] == 2){
        dataIn[assessRow, "Tri"] <- dataIn[assessRow, "Tri"] + 1
      }else if((currRow[1, currMin] == 2 & currRow[1, currMaj] == 2) | ((currRow[1, currMin] + currRow[1, currMaj]) == 4)){
        dataIn[assessRow, "Tet"] <- dataIn[assessRow, "Tet"] + 1
      }
    }
    
    #add seg to plot table, depending on size
    for(currType in c("Tet", "Tri", "LOH", "Loss")){
      tempCol <- currRowPos[currType]
      if(dataIn[assessRow, currType] != 0){
        if(dataIn[assessRow, "nloci"] > sizeCutoff){
          if(dataIn[assessRow, "phyloLoc"] == "T" | dataIn[assessRow, "phyloLoc"] == "BiT"){
            #CNV is large clonal
            plotLargeGains[largeCNVCounter[paste("Clonal", currType, sep="")], tempCol] <- dataIn[assessRow, "nloci"]
            largeCNVCounter[paste("Clonal", currType, sep="")] <- largeCNVCounter[paste("Clonal", currType, sep="")] + 1
          }else if(dataIn[assessRow, "phyloLoc"] == "B"){
            #CNV is large branch
            plotLargeGains[largeCNVCounter[paste("B", currType, sep="")], (tempCol+1)] <- dataIn[assessRow, "nloci"]
            largeCNVCounter[paste("B", currType, sep="")] <- largeCNVCounter[paste("B", currType, sep="")] + 1
          }else{
            #CNV is large leaf
            plotLargeGains[largeCNVCounter[paste("L", currType, sep="")], (tempCol+2)] <- dataIn[assessRow, "nloci"]
            largeCNVCounter[paste("L", currType, sep="")] <- largeCNVCounter[paste("L", currType, sep="")] + 1
          }
        }else{
          if(dataIn[assessRow, "phyloLoc"] == "T"){
            #CNV is large clonal
            plotSmallGains[smallCNVCounter[paste("Clonal", currType, sep="")], tempCol] <- dataIn[assessRow, "nloci"]
            smallCNVCounter[paste("Clonal", currType, sep="")] <- smallCNVCounter[paste("Clonal", currType, sep="")] + 1
          }else if(dataIn[assessRow, "phyloLoc"] == "B"){
            #CNV is large branch
            plotSmallGains[smallCNVCounter[paste("B", currType, sep="")], (tempCol+1)] <- dataIn[assessRow, "nloci"]
            smallCNVCounter[paste("B", currType, sep="")] <- smallCNVCounter[paste("B", currType, sep="")] + 1
          }else{
            #CNV is large leaf
            plotSmallGains[smallCNVCounter[paste("L", currType, sep="")], (tempCol+2)] <- dataIn[assessRow, "nloci"]
            smallCNVCounter[paste("L", currType, sep="")] <- smallCNVCounter[paste("L", currType, sep="")] + 1
          }
        }
      }else{
        #current segmentation doesn't have this type of CNV
        next
      } 
    }
  }
}

############ 2.perform stats ###############
#stats to compare large CNVs
stat1 <- ks.test(plotLargeGains[["ClonalTet"]], plotLargeGains[["BTet"]])$p.value
stat2 <- ks.test(plotLargeGains[["ClonalTet"]], plotLargeGains[["LTet"]])$p.value
stat3 <- ks.test(plotLargeGains[["ClonalTri"]], plotLargeGains[["BTri"]])$p.value
stat4 <- ks.test(plotLargeGains[["ClonalTri"]], plotLargeGains[["LTri"]])$p.value
stat5 <- ks.test(plotLargeGains[["ClonalLOH"]], plotLargeGains[["BLOH"]])$p.value
stat6 <- ks.test(plotLargeGains[["ClonalLOH"]], plotLargeGains[["LLOH"]])$p.value
stat7 <- ks.test(plotLargeGains[["ClonalLoss"]], plotLargeGains[["BLoss"]])$p.value
stat8 <- ks.test(plotLargeGains[["ClonalLoss"]], plotLargeGains[["LLoss"]])$p.value

#stats to compare small CNVs
stat1s <- ks.test(plotSmallGains[["ClonalTet"]], plotSmallGains[["BTet"]])$p.value
stat2s <- ks.test(plotSmallGains[["ClonalTet"]], plotSmallGains[["LTet"]])$p.value
stat3s <- ks.test(plotSmallGains[["ClonalTri"]], plotSmallGains[["BTri"]])$p.value
stat4s <- ks.test(plotSmallGains[["ClonalTri"]], plotSmallGains[["LTri"]])$p.value
stat5s <- ks.test(plotSmallGains[["ClonalLOH"]], plotSmallGains[["BLOH"]])$p.value
stat6s <- ks.test(plotSmallGains[["ClonalLOH"]], plotSmallGains[["LLOH"]])$p.value
stat7s <- ks.test(plotSmallGains[["ClonalLoss"]], plotSmallGains[["BLoss"]])$p.value
stat8s <- ks.test(plotSmallGains[["ClonalLoss"]], plotSmallGains[["LLoss"]])$p.value


############ 3.plot graphs clonal vs subclonal ###############
plotCols <- c(rep("red4", 3), rep("red", 3), rep("blue", 3), rep("lightsteelblue", 3))
pdf(file=paste(subSample[1,6], CNVout, "CNVs/CRC.segSizeDist.pdf", sep=""), width=10, height=8)
  par(xpd=TRUE, mar=c(8,5,5,5), mfrow=c(2,1))
  
  #### large CNVs ####
  boxplot(plotLargeGains, las=2, xaxt='n', pch=20, main=paste("CNV size distributions (> ", sizeCutoff, "Mb)", sep=""), col=plotCols)
  axis(side=1, at=c(1:12), labels = c("ClonalTet", "BTet", "LTet", "ClonalTri", "BTri", "LTri", "ClonalLOH", "BLOH", "LLOH", "ClonalLoss", "BLoss", "LLoss"), las=2)
  
  #add number of CNVs in each distribution
  for(currLen in 1:12){
    nLab <- length(plotLargeGains[!is.na(plotLargeGains[currLen]), 1])
    text(currLen, 275, paste("n = ", nLab, sep=""), cex=0.8)
  }
  
  #add p-values to graph 
  text(2, 295, paste("p = ", signif(stat1, digits = 3)), cex=0.5)
  text(3, 295, paste("p = ", signif(stat2, digits = 3)), cex=0.5)
  text(5, 295, paste("p = ", signif(stat3, digits = 3)), cex=0.5)
  text(6, 295, paste("p = ", signif(stat4, digits = 3)), cex=0.5)
  text(8, 295, paste("p = ", signif(stat5, digits = 3)), cex=0.5)
  text(9, 295, paste("p = ", signif(stat6, digits = 3)), cex=0.5)
  text(11, 295, paste("p = ", signif(stat7, digits = 3)), cex=0.5)
  text(12, 295, paste("p = ", signif(stat8, digits = 3)), cex=0.5)
  
  
  #### small CNVs ####
  boxplot(plotSmallGains, las=2, xaxt='n', pch=20, main=paste("CNV size distributions (< ", sizeCutoff, "Mb)", sep=""), col=plotCols)
  axis(side=1, at=c(1:12), labels = c("ClonalTet", "BTet", "LTet", "ClonalTri", "BTri", "LTri", "ClonalLOH", "BLOH", "LLOH", "ClonalLoss", "BLoss", "LLoss"), las=2)
  
  #add number of CNVs in each distribution
  for(currLen in 1:12){
    nLab <- length(plotSmallGains[!is.na(plotSmallGains[currLen]), 1])
    text(currLen, (sizeCutoff*1.1), paste("n = ", nLab, sep=""), cex=0.8)
  }
  
  #add p-values to graph 
  text(2, (sizeCutoff*1.2), paste("p = ", signif(stat1s, digits = 3)), cex=0.5)
  text(3, (sizeCutoff*1.2), paste("p = ", signif(stat2s, digits = 3)), cex=0.5)
  text(5, (sizeCutoff*1.2), paste("p = ", signif(stat3s, digits = 3)), cex=0.5)
  text(6, (sizeCutoff*1.2), paste("p = ", signif(stat4s, digits = 3)), cex=0.5)
  text(8, (sizeCutoff*1.2), paste("p = ", signif(stat5s, digits = 3)), cex=0.5)
  text(9, (sizeCutoff*1.2), paste("p = ", signif(stat6s, digits = 3)), cex=0.5)
  text(11, (sizeCutoff*1.2), paste("p = ", signif(stat7s, digits = 3)), cex=0.5)
  text(12, (sizeCutoff*1.2), paste("p = ", signif(stat8s, digits = 3)), cex=0.5)
dev.off()







############ 3.get modelled CNVs from parsimony model ###############

parseModelDir <- "8.CNVphylogenetics/CNVphylo/"
modelNames <- ".phyloCNVs.anno.csv"

modCNVtab <- data.frame(matrix(NA, ncol=6, nrow=1000))
names(modCNVtab) <- c("RegionalGain", "ClonalGain", "RegionalLoss", "ClonalLoss", "RegionalLOH", "ClonalLOH") 

#row locations and counters for plotLargeGains and plotSmallGains tables 
rowLocs <- c(1,3,5)
names(rowLocs) <- c("Gain", "LOH", "BLOH")
modCNVCounter <- c(1,1,1,1,1,1)

for(currSam in 1:length(setNames)){
  
  subSample <- subset(sampleList, sampleList[1]==setNames[currSam])
  sampleNames <- subSample[[2]]
  sampleNames <- sampleNames[-(subSample[1,7]+1)]
  
  #get file
  modFileIn <- read.csv(file=paste(subSample[1,6], parseModelDir, subSample[1,1], modelNames, sep=""), header=FALSE)
  names(modFileIn) <- c("ID", "chrom", "pos", "nloci", "pos2", "minor", "major",  sampleNames, "phyloLoc", "model", "genes", "loci")
  
  #get interval sizes
  modFileIn[4] <- modFileIn[[5]] - modFileIn[[3]]
  
  #assess each row and populate table
  currType <- 0
  for(assessRow in 1:nrow(modFileIn)){
    #assess mutation type
    if(modFileIn[assessRow, "minor"] == 0 & modFileIn[assessRow, "major"] == 1){
      #this is a chromosome loss
      currType <- rowLocs["LOH"]
    }else if(modFileIn[assessRow, "minor"] == 0 & modFileIn[assessRow, "major"] == 2){
      #this is an LOH
      currType <- rowLocs["BLOH"]
    }else{
      #this is a chromosome gain
      currType <- rowLocs["Gain"]
    }
    
    if(modFileIn[assessRow, "phyloLoc"] == "T"){
      currType <- as.numeric(currType) + 1
      modCNVtab[modCNVCounter[currType], currType] <- modFileIn[assessRow, "nloci"]
      modCNVCounter[currType] <- modCNVCounter[currType] + 1
    }else{
      modCNVtab[modCNVCounter[currType], currType] <- modFileIn[assessRow, "nloci"]
      modCNVCounter[currType] <- modCNVCounter[currType] + 1
    }
  }  
}

#rearrange columns
modCNVtab <- modCNVtab[c(2,1,4,3,6,5)]
modCNVtab <- modCNVtab / 1000000


############ 4.plot modelled CNV sizes ###############
pdf(file=paste(subSample[1,6], CNVout, "gainLossPlots/modelled.segSizeDist.pdf", sep=""), width=5, height=5)
par(xpd=TRUE, mar=c(8,5,5,5))

#plot large CNVs

#stats to compare gains against losses
stat1 <- wilcox.test(modCNVtab[["RegionalGain"]], modCNVtab[["RegionalLoss"]])
stat2 <- wilcox.test(modCNVtab[["RegionalGain"]], modCNVtab[["RegionalLOH"]])

#compare clonal-subclonal
stat3 <- wilcox.test(modCNVtab[["RegionalGain"]], modCNVtab[["ClonalGain"]])
stat4 <- wilcox.test(modCNVtab[["RegionalLoss"]], modCNVtab[["ClonalLoss"]])
stat5 <- wilcox.test(modCNVtab[["RegionalLOH"]], modCNVtab[["ClonalLOH"]])

boxplot(modCNVtab, las=2, xaxt='n', pch=20, main="aberration size distributions", outcol="grey85")
axis(side=1, at=c(1:6), labels = c("clonal-gain", "subclonal-gain", "clonal-loss", "subclonal-loss", "clonal-BLOH", "subclonal-BLOH"), las=2)
text(6, 240, paste("gain clonal vs subclonal", signif(stat3$p.value, digits = 2)), cex=0.5)
text(6, 230, paste("loss clonal vs subclonal", signif(stat4$p.value, digits = 2)), cex=0.5)
text(6, 220, paste("BLOH clonal vs subclonal", signif(stat5$p.value, digits = 2)), cex=0.5)

text(6, 210, paste("gain size vs loss size", signif(stat1$p.value, digits = 2)), cex=0.5)
text(6, 200, paste("gain size vs LOH size", signif(stat2$p.value, digits = 2)), cex=0.5)

dev.off()


