# 5.2.sampleROIcoverage

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
sampleList <- read.csv(file="~/PhD/CRCproject/masterLynchList.indels.csv", header=FALSE, stringsAsFactors=FALSE)


#prepended names
prepName <- "_statistics.sample_interval_summary"
holdingDir <- "3.samplingAnalysis/coverageFiles/"
apocritaDir <- "processedBams/subSample/sampledRegions"

#remove unwanted samples
sampleList <- sampleList[sampleList[[1]]!="Polyp.08", ]

setNames <- unique(sampleList[[1]])

#make new summary table
regionsTable <- data.frame(matrix(NA, nrow=230418, ncol=nrow(sampleList)))
names(regionsTable) <- sampleList[[2]]

noSamples <- nrow(sampleList)

#main loop to retrieve depths and mean coverage
for(currSet in 1:noSamples){
  repName <- sampleList[currSet, 2]
  fileInName <- paste(sampleList[currSet, 6], holdingDir, sampleList[currSet, 1], "/", repName, prepName, sep="")
  covIn <- read.table(file=fileInName, sep="\t", stringsAsFactors = FALSE, header=TRUE)
  regionsTable[currSet] <- covIn[["average_coverage"]]
}

#name rows by loci
row.names(regionsTable) <- covIn[["Target"]]

#check for interval errors
# errorRow <- c()
# errCounter <- 1
# for(errorCheck in 1:nrow(regionsTable)){
#   intervalCheck <- strsplit(rownames(regionsTable)[errorCheck], "-")
#   if(length(intervalCheck[[1]])!=2){
#     errorRow[errCounter] <- errorCheck
#     errCounter <- errCounter + 1
#   }
# }
# regionsTable <- regionsTable[-errorRow, ]
# readsTable <- readsTable[-errorRow]


###### 2.assess each region for coverage consistency > 9X ######
keepList <- c()
counter <- 1
for(currRow in 1:nrow(regionsTable)){
  assessRow <- table(regionsTable[currRow, ] >= 9)
  
  #if all samples are covered to greater than 20X, keep this region
  if("TRUE" %in% names(assessRow)){
    if(as.integer(assessRow["TRUE"])==noSamples){
      keepList[counter] <- currRow
      counter <- counter + 1
    }
  }else{
    #discard (do nothing)
    next
  }
}

#subset regions table
subRegionsTable <- regionsTable[keepList, ]
subRegionsName <- paste(sampleList[currSet, 6], holdingDir, "allSamples.subRegions.txt", sep="")
write.table(subRegionsTable, file=subRegionsName, sep="\t", row.names=TRUE, col.names = TRUE, quote = FALSE)

#make new regions file (splitting data to columns as neccesary)
RIOtab <- as.data.frame(x=row.names(subRegionsTable), ncol=1, nrow=nrow(subRegionsTable))
platypusTab <- RIOtab
columnSplit <- data.frame(do.call(rbind, strsplit(as.vector(RIOtab[[1]]), split = ":")))
columnSplit2 <- data.frame(do.call(rbind, strsplit(as.vector(columnSplit[[2]]), split = "-")))
RIOtab[1] <- columnSplit[[1]]
RIOtab[2] <- columnSplit2[[1]]
RIOtab[3] <- columnSplit2[[2]]

newRIOname <- paste(sampleList[currSet, 6], holdingDir, "SeqCap.EZ.subsetted.bed", sep="")
write.table(RIOtab, file=newRIOname, sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)

#calling platypus regions file
platypusFile <- paste("/Users/cross01/PhD/CRCproject/", holdingDir, "/SeqCap.EZ.subsetted.platypus.txt", sep="")
write.table(platypusTab, file=platypusFile, sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)



###### 3. subset regions to separate files ######

holdingDir <- "3.samplingAnalysis/sampledRegions/"
shellDir <- "3.samplingAnalysis/platypusScripts/"
apoBamLoc <- "3.samplingAnalysis/"

#take 10000 samples of these ROI, 100 times
for(currSam in 1:100){ 
#   sampleNum <- sample(c(1:nrow(RIOtab)), 10000)
#   sampleNum <- sampleNum[order(sampleNum)]
#   
#   RIOtabSam <- RIOtab[sampleNum, ]
#   platypusTabSam <- platypusTab[sampleNum, ]
#   
#   #output new files
#   newRIOname <- paste(sampleList[1, 6], holdingDir, "SeqCap.EZ.subsetted.", currSam,".bed", sep="")
#   write.table(RIOtabSam, file=newRIOname, sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
#   
#   platypusFile <- paste(sampleList[1, 6], holdingDir, "SeqCap.EZ.subsetted.platypus.", currSam,".txt", sep="")
#   write.table(platypusTabSam, file=platypusFile, sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)
#   
   apoRIOFile <- paste(sampleList[1, 5], holdingDir, "SeqCap.EZ.subsetted.platypus.", currSam,".txt", sep="")
#   
  #make platypus .sh scripts for each sample set
  for(currSet in 1:length(setNames)){
    currSamName <- setNames[currSet]
    platypusOutFile <- paste(sampleList[currSet, 6], shellDir, "runPlatypus.", currSamName, ".", currSam,".sh", sep="")
    vcfOutName <- paste(sampleList[currSet, 5], "1.platypusCalls/samplingAnalysis/", currSamName, "/", currSamName, ".", currSam,".merged.vcf", sep="")
    logFileName <- paste(sampleList[currSet, 5], "1.platypusCalls/samplingAnalysis/", currSamName, "/", currSamName, ".", currSam,".log", sep="")
    
    totalStringsSH <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 4            # Request 48 CPU cores
#$ -l h_rt=48:0:0      # Request 48 hour runtime
#$ -l h_vmem=10G    # Request 3GB RAM / core, i.e. <24GB total
  
#call variants using platypus (by chromosome)
python Software/Platypus_0.5.2/Platypus.py callVariants \\
--bamFiles=", sampleList[1,5], apoBamLoc, currSamName, ".subset.bamList.txt \\
--regions=", apoRIOFile," \\
--output=", vcfOutName," \\
--refFile=/data/BCI-EvoCa/william/referenceHG19/hs37d5.fa \\
--nCPU=4 \\
--maxVariants=75 \\
--mergeClusteredVariants=1 \\
--minMapQual=1 \\
--bufferSize=25 \\
--logFileName=", logFileName 
, sep="")
    
    lapply(totalStringsSH, write, platypusOutFile, append=FALSE)
  }
  
}


