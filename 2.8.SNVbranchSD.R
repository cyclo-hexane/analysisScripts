# 

######################## notes ########################
# assesses branch shape using the SD of branch length as a % of mean branch length (SD-BL)
#
#
######################## libraries ########################
library(e1071)   
library(moments)

######################## subroutines ########################


######################## main program ########################


#input sampleList from commandline arguments
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- "1.platypusCalls/somaticTotal.0.01/"
namePrepended <- "branch&leaf.txt"
rootPrep <- "root.txt"
outDir <- "2.phylogenetics/"

sampleNames <- unique(sampleList[[1]])
#sampleNames <- sampleNames[-c(20:23)]

#output table
metTable <- data.frame(matrix(NA, nrow=length(sampleNames), ncol=9))
metTable[1] <- sampleNames
names(metTable) <- c("names", "noBiopsies", "skewness", "agostino", "trunkLength", "meanBranchLength", "meanOfTrunk", "SD", "SD-BL")

#calculate SD-BL
for(j in 1:length(sampleNames)){
	print(paste("#### making plot for sample ", sampleNames[j], " ####",sep=""))
	
	#subset main list
	subSample <- subset(sampleList, sampleList[1]==sampleNames[j])
	
	setName <- unique(subSample[[1]])
	noSamples <- subSample[1,8]
	normalIndex <- 1+subSample[1,7]
	samNames <- subSample[[2]]
	normalName <- samNames[normalIndex]
	samNamesNoNorm <- samNames[-normalIndex]

	metTable[j, "noBiopsies"] <- length(samNamesNoNorm)
	if(length(samNamesNoNorm) < 3){
	  next()
	}
	
	#get root length
	rootFileName <- paste(subSample[1,6], holdingDir, subSample[1,1], "/", subSample[1,1], ".", rootPrep, sep="")
	rootData <- read.table(file=rootFileName, sep="\t", header = FALSE)
	rootLength <- nrow(rootData)
	metTable[j, "trunkLength"] <- rootLength
	
	#get mean and SD
	branchLengthList <- c()
	for(currBio in 1:length(samNamesNoNorm)){
	  tempFileName <- paste(subSample[1,6], holdingDir, subSample[1,1], "/", samNamesNoNorm[currBio], ".", namePrepended, sep="")
	  tempData <- read.table(file=tempFileName, sep="\t", header = FALSE)
	  branchLengthList[currBio] <- nrow(tempData)
	}
	metTable[j, "skewness"] <- skewness(x=branchLengthList)
	if(length(branchLengthList) > 7){
	  metTable[j, "agostino"] <- agostino.test(x = branchLengthList)$p.value
	}else{
	  metTable[j, "agostino"] <- NA
	}
	metTable[j, "meanBranchLength"] <- mean(branchLengthList)
	metTable[j, "meanOfTrunk"] <- mean(branchLengthList) / rootLength
	metTable[j, "SD"] <- sd(branchLengthList)
	metTable[j, "SD-BL"] <- sd(branchLengthList) / mean(branchLengthList)
		
}	

#write table
outFile <- paste(subSample[1,6], outDir, "branchStandardDev.txt", sep="")
write.table(metTable, file=outFile, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

