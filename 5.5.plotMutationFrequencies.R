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
outDir <- "3.samplingAnalysis/mutationFrequencies/"

setNames <- unique(sampleList[[1]])
setNames <- setNames[-c(17:20, 24)]
orderList <- c(2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,3,3,3)

#sampled interval sizes (from script 5.3.processSampledVCFs)
intFile <- paste(sampleList[1,6], regionsDir,"SeqCap.EZ.totalSizeOfIntervals.txt", sep="")
intervalSizes <- read.table(file=intFile, sep="\t", header=TRUE)

mutFreqMainTab <- data.frame(matrix(NA, ncol=length(setNames), nrow = 20))
names(mutFreqMainTab) <- setNames

#filter each sampled vcf
for(j in 1:length(setNames)){
  subSample <- subset(sampleList, sampleList[1]==setNames[j])
  sampleNames <- subSample[[2]]
  normalName <- sampleNames[subSample[1,7]+1]
  noSamples <- subSample[1,8]
  sampleNamesNoNor <- sampleNames[sampleNames != normalName]
  
  #diversity table (uses sampled interval sizes to get diversity per Mb)
  mutRateTab <- data.frame(matrix(NA, ncol=length(sampleNamesNoNor), nrow=100))
  names(mutRateTab) <- sampleNamesNoNor
  
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
    
    #get mutation frequency for each sample
    for(currSam in 1:length(sampleNamesNoNor)){
      tempTab <- table(confFile[[paste(sampleNamesNoNor[currSam], ".NV", sep="")]] > 0)
      if("TRUE" %in% names(tempTab)){
        mutRateTab[currIter, sampleNamesNoNor[currSam]] <- as.numeric(tempTab["TRUE"]) / intervalSizes[currIter, 1]
      }
    }
  }
  
  #get mean values
  meanValues <- c()
  for(getMean in 1:ncol(mutRateTab)){
    meanValues[getMean] <- median(mutRateTab[[getMean]], na.rm = TRUE)
  }
  
  #order mutations
  mutRateTab <- mutRateTab[order(meanValues)]
  
  #store mean mutation frequency
  mutFreqTab <- c()
  for(getMean in 1:ncol(mutRateTab)){
    mutFreqTab[getMean] <- mean(mutRateTab[[getMean]], na.rm = TRUE)
  }
  mutFreqMainTab[c(1:length(mutFreqTab)), j] <- mutFreqTab
   
    
  #plot graph
  pdf(file = paste(sampleList[1,6], outDir, setNames[j], ".mutFreq.pdf", sep=""), width = (noSamples/2), height = 5)
    boxplot(mutRateTab, las=2)
  dev.off()
  
}


#plot barchart of mean mutation frequencies
meanList <- c()
for(k in 1:ncol(mutFreqMainTab)){
  meanList[k] <- mean(mutFreqMainTab[[k]], na.rm = TRUE)
}


mutFreqMainTab <- mutFreqMainTab[order(orderList, meanList)]
pdf(file = paste(sampleList[1,6], outDir, "adenomaCarcinoma.mutFreq.pdf", sep=""), width = 5, height = 6)
  boxplot(mutFreqMainTab[1:15], las=2)
dev.off()

pdf(file = paste(sampleList[1,6], outDir, "MSIlynch.mutFreq.pdf", sep=""), width = 5, height = 8)
boxplot(mutFreqMainTab[16:19], las=2)
dev.off()

