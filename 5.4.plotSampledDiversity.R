# 5.4.plotSampleDiversity

################### notes: #####################

# processing order:
#
#
# 5.0.countBamReads [run on Apocrita, for each set] (gets tables of read numbers for a given bed file )
#        |
#        v
# 5.1.subsetBamFiles (calculates proportion to be subsampled for each set and makes .sh scripts to subset bams) 
#        |
#        v
# coverageBed run on samples for specified RIO as part of above script
#        |
#        v
# 5.2.sampleROIcoverage (assessed coverage then samples consistent regions 100 times for 10,000 regions)
#        |
#        v
# runPlatypus on sampled regions for subsetted bams
#        |
#        v
# 5.3.processSampledVCFs (to get somatic variants and count interval sampled by each iteration)
#        |
#        v
# 5.4.plotSampleDiversity
#
#################### libraries #################


#################### main program #################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
holdingDir <- "3.samplingAnalysis/platypusCalls/"
regionsDir <- "3.samplingAnalysis/sampledRegions/"
vcfName <- ".somatic.txt"

setNames <- unique(sampleList[[1]])
setNames <- setNames[-c(17:20)]
orderList <- c(2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,3,3,3,3)

#diversity table (uses sampled interval sizes to get diversity per Mb)
diversityTab <- data.frame(matrix(NA, ncol=length(setNames), nrow=100))
names(diversityTab) <- setNames

#sampled interval sizes (from script 5.3.processSampledVCFs)
intFile <- paste(sampleList[1,6], regionsDir,"SeqCap.EZ.totalSizeOfIntervals.txt", sep="")
intervalSizes <- read.table(file=intFile, sep="\t", header=TRUE)

#filter each sampled vcf
for(j in 1:length(setNames)){
  subSample <- subset(sampleList, sampleList[1]==setNames[j])
  sampleNames <- subSample[[2]]
  normalName <- sampleNames[subSample[1,7]+1]
  noSamples <- subSample[1,8]
  nonNormalCols <- c(1:noSamples)
  nonNormalCols <- nonNormalCols[-(subSample[1,7]+1)]
  
  for(currIter in 1:100){
    print(paste("#### assessing sample ", setNames[j], " iteration ", currIter, " ####",sep=""))
    noSamples <- subSample[1,8]
    confFileName <- paste(subSample[1,6], holdingDir, subSample[1,1],"/", subSample[1,1], ".", currIter, vcfName, sep="")
    if(!file.exists(confFileName)){
      print(paste("iteration ", currIter," was missed"))
      next
    }
    
    confFile <- read.table(file=confFileName, sep="\t", header=FALSE, stringsAsFactors = FALSE)
    names(confFile) <- c("chrom", "pos", "ref", "alt", paste(sampleNames, ".NR", sep=""), paste(sampleNames, ".NV", sep=""))
    
    #remove unwanted samples and normal columns, then correct for blank rows
    remNV <- which(paste(normalName, ".NV", sep="") == names(confFile))
    remNR <- which(paste(normalName, ".NR", sep="") == names(confFile))
    confFile <- confFile[-c(remNR, remNV)]
    if(noSamples > 2){
      #select two samples at random
      samSamples <- sample(nonNormalCols, 2)
      samSamples <- samSamples[order(samSamples)]
      
      #sample names selected
      selectedNames <- sampleNames[samSamples]
      selectedNV <- paste(selectedNames, ".NV", sep="")
      selectedNR <- paste(selectedNames, ".NR", sep="")
      
      subCols <- c("chrom", "pos", "ref", "alt", selectedNR, selectedNV)
      
      confFile <- confFile[subCols]
      
      #remove non-present variants from sampled pair
      remRows <- c()
      remCounter <- 1
      for(currRow in 1:nrow(confFile)){
        assessRow <- confFile[currRow, selectedNV]
        if(sum(assessRow) == 0){
          remRows[remCounter] <- currRow
          remCounter <- remCounter + 1
        }
      }
      if(!is.null(remRows)){
        confFile <- confFile[-remRows, ]
      }
    }else{
      selectedNV <- paste(sampleNames[nonNormalCols], ".NV", sep="") 
    }
    noSamples <- 2
    
    #assess diversity as a measure of variants different per Mb sampled
    divTable <- confFile[selectedNV]
    divTable[divTable > 0] <- 1
    currMbInterval <- intervalSizes[currIter, 1]
    divCounter <- 0
    for(currMut in 1:nrow(divTable)){
      if(sum(divTable[currMut, ]) != 2 & sum(divTable[currMut, ]) != 0){
        divCounter <- divCounter + 1
      }
    }
    
    diversityTab[currIter, j] <- divCounter / currMbInterval 
  }
}

#order samples
meanList <- c()
for(currMean in 1:ncol(diversityTab)){
meanList[currMean] <- mean(diversityTab[[currMean]], na.rm = TRUE)
}
diversityTab <- diversityTab[order(orderList, meanList)]

#remove WGS samples
#diversityTab <- diversityTab[-c(Polyp.08.WGS, Set.10, Set.02, Set.09.Proximal, Set.03, Set.04, Set.09.Distal)]

#plot colours
plotCols <- c(rep("purple", 11), rep("red", 5), rep("steelblue", 4))
plotCols <- plotCols[order(orderList, meanList)]

#produce stats by compiling cancer and adenoma diversity values to two tables
adenomaList <- c("Polyp.02", "Polyp.03", "Polyp.05", "Polyp.08.WGS", "Polyp.09")
#adenomaList <- c("Polyp.02", "Polyp.03", "Polyp.05", "Polyp.09")
adenomaDiv <- c()
startAde <- 1
endAde <- 100
for(currAde in 1:length(adenomaList)){
  adenomaDiv[c(startAde:endAde)] <- diversityTab[[adenomaList[currAde]]]
  startAde <- length(adenomaDiv) + 1
  endAde <- length(adenomaDiv) + length(diversityTab[[adenomaList[currAde]]])
}

cancerList <- c("Set.01", "Set.02", "Set.03", "Set.05", "Set.06", "Set.07", "Set.08", "Set.09.Distal", "Set.09.Proximal", "Set.10")
#cancerList <- c("Set.01", "Set.06", "Set.07", "Set.08")
cancerDiv <- c()
startAde <- 1
endAde <- 100
for(currCan in 1:length(cancerList)){
  cancerDiv[c(startAde:endAde)] <- diversityTab[[cancerList[currCan]]]
  startAde <- length(cancerDiv) + 1
  endAde <- length(cancerDiv) + length(diversityTab[[cancerList[currCan]]])
}

#use ks test to compare diversities
ksPvalue <- wilcox.test(cancerDiv[cancerDiv != 0], adenomaDiv[adenomaDiv != 0])$p.value

#compile samples to one table and perform linear model
# cancerTab <- data.frame(matrix(NA, nrow=(nrow(diversityTab)*10), ncol=2))
# cancerTab[1] <- cancerDiv
# cancerTab[2] <- "cancer"
# adenomaTab <- data.frame(matrix(NA, nrow=(nrow(diversityTab)*5), ncol=2))
# adenomaTab[1] <- adenomaDiv
# adenomaTab[2] <- "adenoma"
# totalTab <- rbind(cancerTab, adenomaTab)
# totalTab <- totalTab[totalTab[[1]]!=0,]
# lmResults <- summary(lm(X1 ~ X2, totalTab))
# 
# #get mean diversity values (exclusing MSI for the cancers)
# cancerList <- c("Set.01", "Set.02", "Set.03", "Set.05", "Set.06", "Set.07", "Set.08", "Set.09.Distal", "Set.09.Proximal", "Set.10")
# cancerDiv <- c()
# startAde <- 1
# endAde <- 100
# for(currCan in 1:length(cancerList)){
#   cancerDiv[c(startAde:endAde)] <- diversityTab[[cancerList[currCan]]]
#   startAde <- length(cancerDiv) + 1
#   endAde <- length(cancerDiv) + length(diversityTab[[cancerList[currCan]]])
# }
# 
# meanCancerDiv <- mean(cancerDiv, na.rm = TRUE)
# upperQrt <- quantile(cancerDiv, na.rm = TRUE)["75%"]
# lowerQrt <- quantile(cancerDiv, na.rm = TRUE)["25%"]
# meanAdenomaDiv <- mean(adenomaDiv, na.rm = TRUE)
# upperAdQrt <- quantile(adenomaDiv, na.rm = TRUE)["75%"]
# lowerAdQrt <- quantile(adenomaDiv, na.rm = TRUE)["25%"]

#remove WGS
#diversityTab <- diversityTab[c("Polyp.02", "Polyp.03", "Polyp.05", "Polyp.09", "Set.01", "Set.06", "Set.07", "Set.08")]
#plotCols <- c("red", "red", "red", "red", "purple", "purple", "purple", "purple")

#plot boxplot of diversity
pdf(file=paste(subSample[1,6], "3.samplingAnalysis/diversityPlots/total.diversityPlot.pdf", sep=""), onefile=TRUE, width=10, height=5)
  par(mar=c(8,4,4,4), xpd=TRUE, mfrow=c(1,2))
  boxplot(diversityTab, main="diversity scores, carcinoma vs adenoma", xlab="", ylab="divergent variants / Mb", xaxt="n", col=plotCols)
  axis(1 ,at=c(1:ncol(diversityTab)), labels=names(diversityTab), lwd=0, las=2, cex.axis=1, line=0.5)
  text(x=2, y=9, labels = paste("mean carcinomas: ", mean(cancerDiv[cancerDiv != 0]), sep=""), cex=0.5)
  text(x=2, y=8, labels = paste("mean adenomas: ", mean(adenomaDiv[adenomaDiv != 0]), sep=""), cex=0.5)
  text(x=2, y=7, labels = paste("Wilcox test (p.value ): ", round(ksPvalue, 8), sep=""), cex=0.5)
  
  #add mean lines
  lines(x=c(0,16), y=c(meanCancerDiv, meanCancerDiv), col="purple", lty=2)
  #lines(x=c(0,16), y=c(upperQrt, upperQrt), col="black", lty=2)
  #lines(x=c(0,16), y=c(lowerQrt, lowerQrt), col="black", lty=2)
  lines(x=c(0,16), y=c(meanAdenomaDiv, meanAdenomaDiv), col="red", lty=2)
  #lines(x=c(0,16), y=c(upperAdQrt, upperAdQrt), col="black", lty=2)
  #lines(x=c(0,16), y=c(lowerAdQrt, lowerAdQrt), col="black", lty=2)
  
  #plot comparison
  boxplot(adenomaDiv[adenomaDiv != 0], cancerDiv[cancerDiv != 0], col=c("red", "purple"))
dev.off()

