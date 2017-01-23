# 5.1.subsetBamFiles

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

#sample list
sampleList <- read.csv(file="~/PhD/CRCproject/3.samplingAnalysis/masterSampleList.sampling.csv", header=FALSE, stringsAsFactors=FALSE)
sampleList <- sampleList[-c(103:109), ]

#prepended names
prepName <- ".readCounts.total.txt"
holdingDir <- "3.samplingAnalysis/readCounts/"

setNames <- unique(sampleList[[1]])

#make new summary table
regionsTable <- data.frame(matrix(NA, nrow=0, ncol=3))
names(regionsTable) <- c("names", "reads", "proportion")

#main loop to retrieve counts
for(currSet in 1:length(setNames)){
  repName <- setNames[currSet]
  fileInName <- paste(sampleList[currSet, 6], holdingDir, repName, prepName, sep="")
  covIn <- read.table(file=fileInName, sep="\t", stringsAsFactors = FALSE, header=TRUE)
  
  regionsTable <- rbind(regionsTable, covIn)
}

minReads <- min(regionsTable[[2]])
regionsTable["proportion"] <- minReads / regionsTable[[2]] 


#make shell scripts to subset bam files
for(currSam in 1:nrow(sampleList)){
  currSamName <- sampleList[currSam, 2]
  currSet <- sampleList[currSam, 1]
  subValue <- regionsTable[regionsTable[[1]]==currSamName, 3]
  
  shellFile <- paste(sampleList[1,6], "runScripts/", currSet, ".", currSamName, ".subsetBam.sh", sep="")
  bamFile <- paste(sampleList[1,5], "processedBams/", sampleList[currSam, 1], "/", currSamName, ".mkdub.bam", sep="")
  regionsBam <- paste(sampleList[1,5], "processedBams/", sampleList[currSam, 1], "/", currSamName, ".regions.bam", sep="")
  newBam <- paste(sampleList[1,5], "processedBams/", sampleList[currSam, 1], "/", currSamName, ".subset.bam", sep="")
  coverageFile <- paste(sampleList[1,5], "processedBams/subSampleCoverage/", sampleList[currSam, 1], "/", currSamName, "_statistics", sep="")
  
  totalStringsSH <- paste("#$ -cwd
#$ -V
#$ -pe smp 1            # Request 48 CPU cores
#$ -l h_rt=48:0:0      # Request 48 hour runtime
#$ -l h_vmem=20G    # Request 3GB RAM / core, i.e. <24GB total
  
#no reads before
echo reads in sample ", currSamName," alignment:
samtools view -c -q 1 -L /data/BCI-EvoCa/william/referenceHG19/sureSelectRegions.UNIX.bed", bamFile,"


#subset variants to regions using bedtools
intersectBed -abam ", bamFile," -b /data/BCI-EvoCa/william/referenceHG19/sureSelectRegions.UNIX.bed > ", regionsBam,"


#randomly sample reads to given depth using picard
/usr/java/jdk1.7.0_04/bin/java -Xmx9G -jar ~/bin/DownsampleSam.jar INPUT=", regionsBam," OUTPUT=", newBam," PROBABILITY=", round(subValue, 2),"


#no reads after
echo reads in subsetted sample ", currSamName," alignment:
samtools view -c -q 1 -L /data/BCI-EvoCa/william/referenceHG19/sureSelectRegions.UNIX.bed ", newBam,"
  
#get coverage for these regions
/usr/java/jdk1.7.0_04/bin/java -Xmx10G -jar bin/GenomeAnalysisTK.jar \\
-I",newBam,"\\
-R /data/BCI-EvoCa/william/referenceHG19/hs37d5.fa \\
-T DepthOfCoverage \\
-o", coverageFile,"\\
-L /data/BCI-EvoCa/william/referenceHG19/sureSelectRegions.UNIX.bed \\
--start 1 \\
--stop 1000 \\
--nBins 999  
")
  lapply(totalStringsSH, write, shellFile, append=FALSE)
}

