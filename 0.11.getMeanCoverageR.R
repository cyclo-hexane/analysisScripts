sampleListName <- "~/PhD/CRCproject/masterSampleList.allSamples.filt.csv"
sampleList <- read.csv(file=sampleListName, stringsAsFactors = FALSE, header = FALSE)
setNames <- unique(sampleList[[1]])

dataOut <- data.frame(matrix(NA, ncol=2, nrow=nrow(sampleList)))
names(dataOut) <- c("set", "depth")

#loop thought sets and 
for(curr in 1:nrow(sampleList)){
  
  currSam <- sampleList[curr, 4]
  covIn <- read.table(file=paste(sampleList[1,6], "0.5.bamFileQC/coverageFiles/", currSam, "_statistics.sample_interval_summary", sep=""), sep="\t", stringsAsFactors = FALSE, header = TRUE)
  
  dataOut[curr, 1] <- currSam
  dataOut[curr, 2] <- round(mean(covIn[["average_coverage"]], na.rm = TRUE), digits = 2)
}

write.csv(dataOut, file="~/PhD/CRCproject/0.5.bamFileQC/coverageFiles/summaryTab.1.csv", row.names = FALSE, col.names = TRUE, quote=FALSE)


