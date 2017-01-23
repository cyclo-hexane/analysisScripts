# assesses CNV diversity via segmentations
# plot %genome at diploid against CNV diversity
 
################### notes ###################


################# libraries #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#remSampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.total.exclusions.csv", header=FALSE, stringsAsFactors=FALSE)

CNVholdingDir <- "7.CNVcalls.final/CNVannotation/"
CNVPrepended <- ".CNV.anno.txt"

sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[-c(21:24)]

countsTab <- data.frame(matrix(NA, ncol=length(sampleNames), nrow=3))
names(countsTab) <- sampleNames
row.names(countsTab) <- c("clonal", "subclonal", "propclonal")

#annotate each set each sample 
for(currSet in 1:length(sampleNames)){
	print(paste("#### getting diversity for ", sampleNames[currSet], " ####",sep=""))
	
	#subset main list and get info
  currentSamples <- subset(sampleList, sampleList[1]==sampleNames[currSet])
  
  normalSample <- which(currentSamples[[2]] %in% currentSamples[currentSamples[1,7]+1, 2])
  currentSamples <- currentSamples[-normalSample, ]
  biopsyNames <- currentSamples[[2]]
  majCols <- paste(biopsyNames, "_Major", sep="")
  minCols <- paste(biopsyNames, "_Minor", sep="")
  noSamples <- nrow(currentSamples)
  setName <- currentSamples[1,1]

	#setup input/output names and get files
	CNVfile <- paste(currentSamples[1,6], CNVholdingDir, setName, CNVPrepended, sep="")
	segCNV <- read.table(file=CNVfile, sep="\t", fill=TRUE, stringsAsFactors=FALSE, header=TRUE)
	
	countsTab["clonal", currSet] <- nrow(segCNV[segCNV[["phyloLoc"]]=="T" | segCNV[["phyloLoc"]]=="BiT", ])
	countsTab["subclonal", currSet] <- nrow(segCNV[segCNV[["phyloLoc"]]!="T" & segCNV[["phyloLoc"]]!="BiT", ])
	countsTab["propclonal", currSet] <- countsTab["clonal", currSet]  / nrow(segCNV)
	
}



