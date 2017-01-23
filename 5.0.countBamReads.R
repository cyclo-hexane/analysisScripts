# 5.0.countBamReads [run on APOCRITA]

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


# to run this sript make an Apocrita .sh file for each set containing and R run command with arguments stated below

arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
if(length(arguments)!=2){
  stop("\n#### please use syntax > Rscript 5.0.countBamReads.R <sample list file> <sample set name> ####\n")
}

sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="/data/BCI-EvoCa/william/CRCproject/3.samplingAnalysis/masterBamList.csv", header=FALSE, stringsAsFactors=FALSE)

currentSample <- arguments[2]
setNames <- subset(sampleList, sampleList[1]==currentSample)

holdingDir <- "processedBams/"

#read in coverage reference table
#coverageTab <- read.table(file="/data/BCI-EvoCa/william/CRCproject/3.samplingAnalysis/allSamples.subRegions.txt", header=TRUE, stringsAsFactors=FALSE)
#ROIfile <- read.table(file="/data/BCI-EvoCa/william/CRCproject/3.samplingAnalysis/subRegions.txt", header=FALSE, stringsAsFactors=FALSE)

readsTab <- data.frame(matrix(NA, ncol=3, nrow=nrow(setNames)))
readsTab[1] <- setNames[[2]]
names(readsTab) <- c("names", "reads", "proportion")

for(currSam in 1:nrow(setNames)){
  print(paste("#### ", setNames[currSam, 2], "####"))
  
  #bam name
  bamFile <- paste(setNames[currSam, 5], holdingDir, setNames[currSam, 1], "/", setNames[currSam, 2], ".mkdub.bam", sep="")
  
  #run samtools and count reads in region
  noReads <- system(paste("samtools view -c -q 1 -L /data/BCI-EvoCa/william/referenceHG19/sureSelectRegions.UNIX.bed", bamFile), intern = TRUE)
  
  readsTab[currSam, 2] <- as.numeric(noReads)
}

minReads <- min(readsTab[[2]])
readsTab[3] <- minReads / readsTab[[2]]

#write table
write.table(readsTab, file=paste(setNames[currSam, 5], holdingDir, setNames[1, 1], ".readCounts.total.txt", sep=""), sep="\t", row.names=FALSE, quote = FALSE)

