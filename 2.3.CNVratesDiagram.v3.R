# plot summary of CNV disruptions from cloneHD segmentations

##########################   notes   ##########################
#         1. plot cross-sample summary of CNVs
#         2. calculate cross-cancer heterogeneity
#

########################## libraries  ##########################

######################### subroutines ##########################

######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[[1]])

#loop through samplelist and remove normals
normalCols <- c()
normalCounter <- 1
for(currSet in 1:length(setNames)){
  currTab <- sampleList[sampleList[[1]]==setNames[currSet], ]
  tempNorCol <- which(currTab[(1 + currTab[1,7]), 2] == sampleList[[2]])
  for(currNo in 1:length(tempNorCol)){
    normalCols[normalCounter] <- tempNorCol[currNo] 
    normalCounter <- normalCounter + 1 
  }
}
normalCols <- unique(normalCols)
sampleList <- sampleList[-normalCols, ]

#working dir
CNVholdingDir <- "7.CNVcalls.final/baf.updated/"
CNVfileNames <- ".penalty0.95.baf.gt.txt"

plotDir <- "9.driverEvents/CNVs/"
percentDisr <- c()
sampleNames <- unique(sampleList[[1]])
noSets <- length(sampleNames)

#chromosome positions
chromPos <- read.csv(file="~/PhD/ReferenceGenome/chromosomeStartEnd.csv", stringsAsFactors = FALSE, header=TRUE)
totalGenLength <- sum(as.numeric(chromPos[["size"]]))

#sample position (in the list)
samplePos <- c()

#remove samples from list
#sampleList <- sampleList[1:73,]
#setNames <- setNames[c(1:11)]
#sampleList <- sampleList[c(74:79,83:96,106:113),]
#setNames <- setNames[c(12:16, 20:23)]
sampleList <- sampleList[1:96,]
setNames <- setNames[c(1:16)]

#get annotation table with gene names
cnvTab <- read.csv(file="~/PhD/CRCproject/9.driverEvents/CNVs/CNVfreqTable.csv", stringsAsFactors = FALSE, header = TRUE)
cnvTab <- cnvTab[cnvTab[["gene"]]!="N", ]
cnvTab <- cnvTab[, -c(6:ncol(cnvTab))]
names(cnvTab)[3:5] <- c("chrom", "gPosStart", "gPosEnd")

#table to store clonal\subclonal CNA %
CNAsubTable <- data.frame(matrix(NA, nrow = nrow(sampleList), ncol = 4))
names(CNAsubTable) <- c("gain", "subgain", "loss", "subloss")
row.names(CNAsubTable) <- sampleList[[2]]

#get loci positions and alter to get true genome positions
cytobandList <- read.table(file="~/PhD/ReferenceGenome/cytoBand.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
for(annGene in 1:nrow(cnvTab)){
  splitChrom <- cnvTab[annGene, 1]
  if(length(which(strsplit(splitChrom, "")[[1]]=="p")) != 0){
    lociTemp <- paste("p", strsplit(splitChrom, split = "p")[[1]][2], sep="")
    splitChrom <- strsplit(splitChrom, split = "p")[[1]][1]
    subChromPos <- cytobandList[cytobandList[[1]]==splitChrom, ]
    centromerePos <- subChromPos[subChromPos[[5]]=="acen", ]
    subChromPos <- subChromPos[subChromPos[[2]] <= centromerePos[1, 2], ]
  }else{
    lociTemp <- paste("q", strsplit(splitChrom, split = "q")[[1]][2], sep="")
    splitChrom <- strsplit(splitChrom, split = "q")[[1]][1]
    subChromPos <- cytobandList[cytobandList[[1]]==splitChrom, ]
    centromerePos <- subChromPos[subChromPos[[5]]=="acen", ]
    subChromPos <- subChromPos[subChromPos[[2]] > centromerePos[2, 3], ]
  }
  cnvTab[annGene, c(3:5)] <- subChromPos[subChromPos[[4]]==lociTemp, c(1:3)]
  addPostionValue <- chromPos[chromPos[["chrom"]]==splitChrom, "start"]
  
  cnvTab[annGene, c(4:5)] <- cnvTab[annGene, c(4:5)] + addPostionValue
}




#### loop through all biopsies individually and get CNV data into data list ####
CNVdataList <- as.list(NA)
for(currSeg in 1:nrow(sampleList)){
  print(paste("#### getting segs for sample ", sampleList[currSeg, 2], " ####",sep=""))
  
  setName <- sampleList[currSeg, 1]
  currSamName <- sampleList[currSeg, 2]
  
  #set names for all samples in current set
  currentSamples <- sampleList[sampleList[[1]]==setName, ]
  samNames <- currentSamples[[2]]
  
  #setup input/output names
  dataIn <- read.table(file=paste(sampleList[currSeg, 6], CNVholdingDir, setName, CNVfileNames, sep=""), sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
  dataIn[dataIn == "X"] <- 1
  
  #get current sample columns
  currMin <- which(paste(currSamName, "_Minor", sep="")==names(dataIn))
  currMaj <- which(paste(currSamName, "_Major", sep="")==names(dataIn))
  
  #get non-diploid segmentations and phylogenetic locations
  phyloLoc <- c()
  stateList <- c()
  remRow <- c()
  remCounter <- 1
  colIndexes <- c(currMin, currMaj)
  for(currRem in 1:nrow(dataIn)){
    assessRow <- dataIn[currRem, colIndexes]
    
    #determine if normal seg
    assessNorm <- table(assessRow == 1)
    if("TRUE" %in% names(assessNorm)){
      if(as.numeric(assessNorm["TRUE"])==2){
        remRow[remCounter] <- currRem
        remCounter <- remCounter + 1
        next
      }
    }
    
    #determine phylogenetic location by using all samples in current set
    stateStrings <- c()
    for(currState in 1:length(samNames)){
      stateTemp <- dataIn[currRem, c(paste(samNames[currState], "_Minor", sep=""), paste(samNames[currState], "_Major", sep=""))]
      stateStrings[currState] <- paste(stateTemp[1], ":", stateTemp[2], sep="")
    }
    stateStrings <- table(stateStrings)
    stateStrings <- stateStrings[order(stateStrings)]
    if(length(stateStrings)==1){
      phyloLoc[currRem] <- "T"
      stateList[currRem] <- names(stateStrings)[1]
    }else{
      phyloLoc[currRem] <- "B/L"
      stateStrings <- stateStrings[names(stateStrings)!="1:1"]
      tempStates <- names(stateStrings[which(max(stateStrings)==stateStrings)])
      stateList[currRem] <- tempStates[1]
    }
  }
  
  #assess potential rows (segs) to be removed
  if(!is.null(remRow)){
    phyloLoc <- phyloLoc[-remRow]
    dataIn <- dataIn[-remRow, c(1:4, colIndexes)]
    stateList <- stateList[-remRow]
  }
  dataIn[ncol(dataIn)+1] <- phyloLoc
  dataIn[ncol(dataIn)+1] <- stateList
  names(dataIn)[c((ncol(dataIn)-1):ncol(dataIn))] <- c("phyloLoc", "state")
  
  #save table to
  dataIn["nloci"] <- dataIn[["last.locus"]] - dataIn[["first.locus"]]
  CNVdataList[[currSeg]] <- dataIn[c("chr", "first.locus", "nloci", "last.locus", "phyloLoc", "state")]
  names(CNVdataList)[currSeg] <- currSamName
  
  #calculate percentage disrupted
  percentDisr[currSeg] <- sum(as.numeric(dataIn[["nloci"]])) / totalGenLength
  
  #save sample set position
  samplePos[currSeg] <- which(setName==setNames)
}
print("#### done making segmentation list ####")

#loop through and get mean % genome disruptions per set
tempPercDist <- c()
addCounter <- 1
for(currPerc in 1:length(unique(samplePos))){
  sumPercs <- mean(percentDisr[which(samplePos == currPerc)])
  for(currAdd in 1:length(percentDisr[which(samplePos == currPerc)])){
    tempPercDist[addCounter] <- sumPercs
    addCounter <- addCounter + 1
  }
}
tempPercDist <- round(tempPercDist, digits = 6)

#order by amount of genome disrupted
for(currOrd in 1:length(sampleNames)){
  currSub <- subset(sampleList, sampleList[1]==sampleNames[currOrd])
  ordPerc <- percentDisr[!is.na(match(sampleList[[2]], currSub[[2]]))]
  addValues <- sort(ordPerc, index.return=TRUE, decreasing = TRUE)$ix
  addValues <- addValues / 1000000
  
  #increment percentage values with small number to allow ordering
  tempPercDist[!is.na(match(sampleList[[2]], currSub[[2]]))] <- tempPercDist[!is.na(match(sampleList[[2]], currSub[[2]]))] + addValues 
}


sampleOrder <- order(tempPercDist)
percentDisr <- percentDisr[sampleOrder]



##################### 1.plot summary of CNVs #####################

#get chromosome density segments first

#chromosome intervals table
chrDensity <- data.frame(matrix(0, ncol=length(seq(0, totalGenLength, 10000000)), nrow=5))
chrDensity[1,] <- seq(0, totalGenLength, 10000000)
rownames(chrDensity) <- c("genomePos", "freqGainCl", "freqGain", "freqLossCl", "freqLoss")

#populate plotting table
for(currSam in sampleOrder){
  print(paste("#### adding biopsy ", names(CNVdataList)[[currSam]], " to data ####", sep=""))
  lossgainFlag <- 0
  xcounter <- 0
  
  currSetTemp <- CNVdataList[[currSam]]
  CNAsubTable[currSam, "gain"] <- nrow(currSetTemp[(currSetTemp[["state"]]!="0:1" & currSetTemp[["state"]]!="0:2") & currSetTemp[["phyloLoc"]]=="T", ])  / nrow(currSetTemp[currSetTemp[["state"]]!="0:1" & currSetTemp[["state"]]!="0:2", ])
  CNAsubTable[currSam, "subgain"] <- nrow(currSetTemp[(currSetTemp[["state"]]!="0:1" & currSetTemp[["state"]]!="0:2") & currSetTemp[["phyloLoc"]]=="B/L", ])  / nrow(currSetTemp[currSetTemp[["state"]]!="0:1" & currSetTemp[["state"]]!="0:2", ])
  CNAsubTable[currSam, "loss"] <- nrow(currSetTemp[(currSetTemp[["state"]]=="0:1" | currSetTemp[["state"]]=="0:2") & currSetTemp[["phyloLoc"]]=="T", ])  / nrow(currSetTemp[currSetTemp[["state"]]=="0:1" | currSetTemp[["state"]]=="0:2", ])
  CNAsubTable[currSam, "subloss"] <- nrow(currSetTemp[(currSetTemp[["state"]]=="0:1" | currSetTemp[["state"]]=="0:2") & currSetTemp[["phyloLoc"]]=="B/L", ])  / nrow(currSetTemp[currSetTemp[["state"]]=="0:1" | currSetTemp[["state"]]=="0:2", ])

  
  for(chromo in 1:22){
    currCentro <- chromPos[chromo, "centromere"]
    xcounter <- chromPos[chromo, "start"]
    currSamSegs <- subset(CNVdataList[[currSam]], CNVdataList[[currSam]]["chr"]==chromo)
    if(nrow(currSamSegs)==0){
      next
    }
    #colours are only used as flags for clonal sorting
    for(currRect in 1:nrow(currSamSegs)){
      if(currSamSegs[currRect, "state"]=="0:2" | currSamSegs[currRect, "state"]=="0:1"){
        if(currSamSegs[currRect, "phyloLoc"] == "T"){
          plotCol <- "salmon"
        }else{
          plotCol <- "lightsalmon"
        }
      }else{
        if(currSamSegs[currRect, "phyloLoc"] == "T"){
          plotCol <- "steelblue"
        }else{
          plotCol <- "lightsteelblue"
        }
      }
      
      xTemp <- currSamSegs[currRect, "first.locus"] + xcounter
      xTempEnd <- currSamSegs[currRect, "last.locus"] + xcounter
      
      #record current segmentation to region
      if(plotCol == "salmon" | plotCol == "lightsalmon"){
        lossgainFlag <- 4
      }else{
        lossgainFlag <- 2
      }
      if(currSamSegs[currRect, "nloci"] > 10000000){
        if(plotCol == "lightsalmon" | plotCol == "lightsteelblue"){
          lossgainFlag <- lossgainFlag + 1
        }
        chrDensity[lossgainFlag, as.numeric(chrDensity[1,]) >= xTemp & as.numeric(chrDensity[1,]) <= xTempEnd] <- 1 + chrDensity[lossgainFlag, as.numeric(chrDensity[1,]) >= xTemp & as.numeric(chrDensity[1,]) <= xTempEnd]
      }
      
    }
  }
}

#get total proportions for histogram of clonal v non-clonal
totalGainClonal <- sum(chrDensity[2, ])
totalGainSub <- sum(chrDensity[3, ])
maxGain <- totalGainClonal + totalGainSub

totalLossClonal <- sum(chrDensity[4, ])
totalLossSub <- sum(chrDensity[5, ])
maxLoss <- totalLossClonal + totalLossSub

totalGainClonal <- sum(chrDensity[2, ] * 2.5/maxGain)
totalGainSub <- sum(chrDensity[3, ] * 2.5/maxGain)
totalLossClonal <- sum(chrDensity[4, ] * 2.5/maxLoss)
totalLossSub <- sum(chrDensity[5, ] * 2.5/maxLoss)

#get proportions to fit density to graph
divNumber <- max(chrDensity[2,])
divNumber2 <- max(chrDensity[3,])
divNumber <- divNumber + divNumber2
chrDensity[2,] <- chrDensity[2,] * (2.5/divNumber)
chrDensity[3,] <- chrDensity[3,] * (2.5/divNumber)
chrDensity[4,] <- chrDensity[4,] * (2.5/divNumber)
chrDensity[5,] <- chrDensity[5,] * (2.5/divNumber)




################## main plot graph and chromosome CNV density #####################

xmax <- max(chromPos["start"])
outputLoc <- paste(sampleList[1,6], plotDir, "totalCNVsummary.carcinoma.pdf", sep="")
pdf(file=outputLoc, width=10, height=(noSets/3))
  #percentage plot start-end
  plotStart <- 2900000000
  plotEnd <- 3000000000
  plotDiff <- plotEnd - plotStart 
  
  par(mar=c(8,5,2,5), xpd=TRUE, cex.main = 2.5) 
  plot(1, 1, col="white", axes=F, xlim=c(0, plotEnd), ylim=c(-10, ((nrow(sampleList)+1)/2)), xlab="", ylab="", main="")
  
  ycounter <- 0
  
  #colour vector (order of MSI, adenoma, cancer)
  colVector <- rep("darkslateblue", nrow(sampleList))
  
  #plot segmentations
  for(currSam in (sampleOrder)){
    xcounter <- 0
    lossgainFlag <- 0
    for(chromo in 1:22){
      currCentro <- chromPos[chromo, "centromere"]
      xcounter <- chromPos[chromo, "start"]
      currSamSegs <- subset(CNVdataList[[currSam]], CNVdataList[[currSam]]["chr"]==chromo)
      if(nrow(currSamSegs)==0){
        next
      }
      for(currRect in 1:nrow(currSamSegs)){
        if(currSamSegs[currRect, "state"]=="0:2" | currSamSegs[currRect, "state"]=="0:1"){
          if(currSamSegs[currRect, "phyloLoc"] == "T"){
            plotCol <- "blue"
          }else{
            plotCol <- "lightsteelblue"
          }
        }else{
          if(currSamSegs[currRect, "phyloLoc"] == "T"){
            plotCol <- "red"
          }else{
            plotCol <- "pink"
          }
        }
        
        xTemp <- currSamSegs[currRect, "first.locus"] + xcounter
        xTempEnd <- currSamSegs[currRect, "last.locus"] + xcounter
        rect(xleft=xTemp, xright=xTempEnd, ybottom=(ycounter+0.4), ytop=ycounter+0.6, col=plotCol, border = NA)
      }
    }
    ycounter <- ycounter + 0.5
  }
  
  #add chromosome partitions and names
  for(currLine in 1:22){
    lines(x=c(chromPos[currLine, "start"], chromPos[currLine, "start"]),  y=c(0, (nrow(sampleList)/2)), lty=1)
    text(x=chromPos[currLine, "centromere"], y=ycounter+0.5, labels = currLine, cex=0.6)
  }
   
  #add sample names
  counter <- 1
  nameCounter <- 1
  for(currName in sampleOrder){
    text(pos = 2, x=c(-500, -500), y=c((counter-0.5), (counter-0.5)), labels=names(CNVdataList)[currName], col="black", cex=0.5)
    counter <- counter + 0.5
    nameCounter <- nameCounter + 1
  }
  
  ycounter2 <- -5
  
  #add density graph per chromosome
  for(currChr in 2:ncol(chrDensity)){
    lossTemp <- chrDensity[4, currChr]
    gainTemp <- chrDensity[2, currChr]
    rect(xleft=chrDensity[1, (currChr-1)], xright=chrDensity[1, currChr], ytop=ycounter2+gainTemp, ybottom=ycounter2, col="red", border = NA)
    rect(xleft=chrDensity[1, (currChr-1)], xright=chrDensity[1, currChr], ytop=ycounter2+gainTemp+chrDensity[3, currChr], ybottom=ycounter2+gainTemp, col="pink", border = NA)
    rect(xleft=chrDensity[1, (currChr-1)], xright=chrDensity[1, currChr], ybottom=ycounter2-lossTemp, ytop=ycounter2, col="blue", border = NA)
    rect(xleft=chrDensity[1, (currChr-1)], xright=chrDensity[1, currChr], ybottom=ycounter2-lossTemp-chrDensity[5, currChr], ytop=ycounter2-lossTemp, col="lightsteelblue", border = NA)
  }
  
  #add chromosome partitions and names
  for(currLine in 1:22){
    lines(x=c(chromPos[currLine, "start"], chromPos[currLine, "start"]),  y=c(-2.5, -7.5), lty=2, cex=0.75)
  }
   
  #add density axis
  lines(x=c(chromPos[1, "start"], chromPos[1, "start"]),  y=c(-2.5, -7.5), lty=1)
  text(pos = 2, x=c(-500, -500), y=ycounter2, labels = "diploid", cex=1)
  text(pos = 2, x=c(-500, -500), y=ycounter2+2.5, labels = "gain", cex=1)
  text(pos = 2, x=c(-500, -500), y=ycounter2-2.5, labels = "loss", cex=1)
  
  #add gene annotations
  #for(addAnno in 1:nrow(cnvTab)){
  #  text(x = cnvTab[addAnno, "gPosStart"], y = -9, labels = cnvTab[addAnno, "gene"], cex=0.5, srt=45)
  #  lines(x = c(cnvTab[addAnno, "gPosStart"], cnvTab[addAnno, "gPosStart"]), y = c(-8, -5))
  #}
  
  #histogram of clonal v non-clonal
  #gains
  xPosTempL <- chrDensity[1, (currChr-1)] + 50000000
  xPosTempR <- xPosTempL + 100000000
  rect(xleft=xPosTempL, xright=xPosTempR, ytop=ycounter2+totalGainClonal, ybottom=ycounter2, col="red", border = NA)
  rect(xleft=xPosTempL, xright=xPosTempR, ytop=ycounter2+totalGainClonal+totalGainSub, ybottom=ycounter2+totalGainClonal, col="pink", border = NA)
  #losses
  rect(xleft=xPosTempL, xright=xPosTempR, ybottom=ycounter2-totalLossClonal, ytop=ycounter2, col="blue", border = NA)
  rect(xleft=xPosTempL, xright=xPosTempR, ybottom=ycounter2-totalLossClonal-totalLossSub, ytop=ycounter2-totalLossClonal, col="lightsteelblue", border = NA)
  xPosTempR <- xPosTempR + 50000000
  #add % marks to histogram
  gainPerc <- round((totalGainClonal / 2.5) * 100, digits = 1)
  text(x=(xPosTempR), y=ycounter2+totalGainClonal, labels=gainPerc, cex=0.5)
  lossPerc <- round((totalLossClonal / 2.5) * 100, digits = 1)
  text(x=(xPosTempR), y=ycounter2-totalLossClonal, labels=lossPerc, cex=0.5)
  
  #add percentage genome disrupted
  counter <- 0
  for(currBar in 1:length(percentDisr)){
    rect(xleft=plotStart, xright=(plotStart+plotDiff*percentDisr[currBar]), ybottom=counter+0.4, ytop=counter+0.6, col=colVector[currBar])
    counter <- counter + 0.5
  }
  
  #add percentage disruption axis
  ycounter <- ycounter + 0.25
  lines(x=c(plotStart, plotEnd),  y=c(ycounter, ycounter), lty=1)
  ycounter <- ycounter + 0.25
  text(labels="0",  x=plotStart, y=ycounter, cex=0.5)
  text(labels="100%",  x=plotEnd, y=ycounter, cex=0.5)
  
dev.off()



#stats to compare adenomas and carcinoma clonal and subclonal CNAs

#stat diff of losses
wilcox.test(CNAsubTable[1:73, 2], CNAsubTable[74:96, 2])

#stat diff of losses
wilcox.test(CNAsubTable[1:73, 4], CNAsubTable[74:96, 4])


