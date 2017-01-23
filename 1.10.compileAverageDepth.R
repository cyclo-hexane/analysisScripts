# takes interval_summary files from GATK coverageBam and compiles a table

################### notes: #####################


#################### libraries #################


#################### main program #################

#sample list
sampleList <- read.csv(file="~/PhD/CRCproject/3.samplingAnalysis/masterSampleList.sampling.csv", header=FALSE, stringsAsFactors=FALSE)

#prepended names
prepName <- "_statistics.sample_interval_summary"
holdingDir <- "bamFileQC/coverageFiles/"

#remove polyp.08 and duplicate normal from Set.09
sampleList <- sampleList[-c(73, 103:109), ]

#make new summary table
regionsTable <- data.frame(matrix(NA, nrow=nrow(sampleList), ncol=2))
names(regionsTable) <- c("av.reads", "av.depth")
rownames(regionsTable) <- sampleList[[2]]

#main loop to retrieve counts
for(currSet in 7:nrow(regionsTable)){
  fileInName <- paste(sampleList[currSet, 6], holdingDir, sampleList[currSet, 2], prepName, sep="")
  covIn <- read.table(file=fileInName, sep="\t", stringsAsFactors = FALSE, header=TRUE)
  
  #remove regions with no coverage (Y chromosome in females)
  covIn <- covIn[!is.na(covIn[[3]]), ]
  
  regionsTable[currSet, 1] <- mean(covIn[[2]])
  regionsTable[currSet, 2] <- mean(covIn[[3]]) 
}

#output table
fileOut <- paste(sampleList[currSet, 6], holdingDir, "averageCoverageFile.csv", sep="")
write.table(regionsTable, file=fileOut, sep=",", row.names = TRUE, col.names = TRUE, quote = FALSE)
