#script checks for the existence of driver mutations in SNP files
#drivers are assessed from CRC project (TCGA/COSMIC) list. Makes output table of suspect variants
#also count number of exonic variants

#get sampleList
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterLynchList.polyps.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- sampleList[sampleList[[1]]=="Ha.53", ]

#get driver list
driverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table.05-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)
#driverList <- read.csv(file="~/PhD/CRCproject/9.driverEvents/supp.table.05-lynchDrivers.csv", header=TRUE, stringsAsFactors=FALSE)

setNames <- unique(sampleList[[1]])
holdingDir <- "1.platypusCalls/exome.0.01/"
#appendedName <- ".snv.annoVar.variant_function.txt"
appendedName <- ".snv.annoVar.exonic_variant_function.0.01.txt"

#to count variants
SNVdir <- "1.platypusCalls/exome.0.01/"
appendedSNV <- ".snv.annoVar.exonic_variant_function.0.01.txt"

#to count non-synon and synon
dndsTable <- data.frame(matrix(NA, nrow = length(setNames), ncol=8))
names(dndsTable) <- c("set", "non-syn.clonal", "synon.clonal", "non-syn.sub", "synon.sub", "clonalRatio", "subClonalRatio", "noBiopsies")

fileoutPrep <- ".driverSNVs.txt"
#fileoutPrep <- ".driverIndels.txt"

#typeFile <- "genome"
typeFile <- "exome"

#main loop
for(currSet in 1:length(setNames)){
  subSample <- subset(sampleList, sampleList[1]==setNames[currSet])
  noSamples <- subSample[1,8]
  sampleNames <- subSample[[2]]
  normalIndex <- subSample[1,7] + 1
  nonNorSamples <- sampleNames[-normalIndex]
  
  #SNPs file
  fileName <- paste(subSample[1,6], holdingDir, subSample[1,1], "/", subSample[1,1], appendedName, sep="")
  dataIn <- read.table(file=fileName, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  names(dataIn) <- c("line", "type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(sampleNames, ".NV", sep=""), sampleNames)
  
  if(typeFile == "genome"){
    dataIn <- as.data.frame(append(dataIn, values = NA, after = 0))
  }
  
  #subset data by non-synon variants
  #dataIn <- dataIn[dataIn[2]=="nonsynonymous SNV" | dataIn[2]=="stopgain" | dataIn[2]=="stoploss", ]
  
  #get phyloLoc vector
  clonalLoc <- c()
  counter <- 1
  for(currPhy in 1:nrow(dataIn)){
    tempTab <- table(dataIn[currPhy, nonNorSamples] > 0)
    if("TRUE" %in% names(tempTab)){
      if(as.numeric(tempTab["TRUE"]) == length(nonNorSamples)){
        clonalLoc[counter] <- currPhy
        counter <- counter + 1  
      }
    }
  }
  
  #count non-synony and synon
  clonalVars <- dataIn[clonalLoc,]
  subclonalVars <- dataIn[-clonalLoc,]
  typeClonalTab <- table(clonalVars[[2]])
  typeSubclonalTab <- table(subclonalVars[[2]])
  if("stopgain" %in% names(typeClonalTab)){
    dndsTable[currSet, "non-syn.clonal"] <- sum(typeClonalTab[c("nonsynonymous SNV", "stopgain")])
  }else{
    dndsTable[currSet, "non-syn.clonal"] <- typeClonalTab["nonsynonymous SNV"]
  }
  dndsTable[currSet, "synon.clonal"] <- sum(typeClonalTab["synonymous SNV"])
  if("stopgain" %in% names(typeSubclonalTab)){
    dndsTable[currSet, "non-syn.sub"] <- sum(typeSubclonalTab[c("nonsynonymous SNV", "stopgain")])
  }else{
    dndsTable[currSet, "non-syn.sub"] <- typeSubclonalTab["nonsynonymous SNV"]
  }
  dndsTable[currSet,  "synon.sub"] <- sum(typeSubclonalTab["synonymous SNV"])
  dndsTable[currSet,  "clonalRatio"] <- dndsTable[currSet, "non-syn.clonal"] / dndsTable[currSet, "synon.clonal"]
  dndsTable[currSet,  "subClonalRatio"] <- dndsTable[currSet, "non-syn.sub"] / dndsTable[currSet, "synon.sub"]
  dndsTable[currSet, "set"] <- setNames[currSet]
  dndsTable[currSet, "noBiopsies"] <- length(nonNorSamples) 
  
  #separate gene names
  for(currSpl in 1:nrow(dataIn)){
    dataIn[currSpl, 1] <- strsplit(as.character(dataIn[currSpl, 3]), split = ":")[[1]][1]
  }
  
  #get driver genes
  keepList <- c()
  counter <- 1
  for(currRow in 1:nrow(dataIn)){
    if(dataIn[currRow, 1] %in% driverList[[1]]){
      keepList[counter] <- currRow
      counter <- counter + 1
    }
  }
  dataIn <- dataIn[keepList, ]
  
  #get allele frequencies
  for(currSam in 1:noSamples){
    dataIn[ncol(dataIn) + 1] <- dataIn[8 + noSamples + currSam] / dataIn[8 + currSam]
  }
  
  #name columns
  names(dataIn) <- c("line", "type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(sampleNames, ".NR", sep=""), paste(sampleNames, ".NV", sep=""), paste(sampleNames, ".AF", sep=""))
  
  #sort by allele frequency in normal
  normalCol <- 8 + noSamples + noSamples + normalIndex
  dataIn <- dataIn[order(dataIn[normalCol]), ]
  
  #output driver query list
  fileOut <- paste(subSample[1,6], holdingDir, subSample[1,1], "/", subSample[1,1], fileoutPrep, sep="")
  write.table(dataIn, file=fileOut, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)  
}




#make summary table
varTab <- data.frame(matrix(NA, ncol=2, nrow=nrow(sampleList)))
names(varTab) <- c("sampleName", "noExomicVar")
varTab[1] <- sampleList[[2]]

#count exonic variants
for(currSet in 1:length(setNames)){
  subSample <- subset(sampleList, sampleList[1]==setNames[currSet])
  noSamples <- subSample[1,8]
  sampleNames <- subSample[[2]]
  normalIndex <- subSample[1,7] + 1
  
  #SNPs file
  fileName <- paste(subSample[1,6], SNVdir, subSample[1,1], "/", subSample[1,1], appendedSNV, sep="")
  dataIn <- read.table(file=fileName, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  names(dataIn) <- c("line", "type", "genes", "chrom", "pos", "pos2", "ref", "alt", paste(subSample[[2]], ".NR", sep=""), subSample[[2]])
  
  sampleNames <- sampleNames[-(subSample[1,7]+1)]
  
  for(currSam in 1:length(sampleNames)){
    tempData <- dataIn[[sampleNames[currSam]]]
    varTab[varTab[[1]]==sampleNames[currSam], 2] <- length(tempData[tempData>0]) 
  }
  
}

#remove normal samples
varTab <- varTab[!is.na(varTab[[2]]), ]

#output driver query list
fileOut <- paste(subSample[1,6], SNVdir, "total.exomeVars.txt", sep="")
write.table(varTab, file=fileOut, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)  



#perform dNdS analysis
write.table(dndsTable, file="~/PhD/CRCproject/1.platypusCalls/exome.0.01/dNdsTab.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)

cancerRatios <- dndsTable[c(1:3,5:11), c("clonalRatio", "subClonalRatio")]
adenomaRatios <- dndsTable[12:20, c("clonalRatio", "subClonalRatio")]

logCan <- log2(abs(cancerRatios[["clonalRatio"]]-cancerRatios[["subClonalRatio"]]))
logAden <- log2(abs(adenomaRatios[["clonalRatio"]]-adenomaRatios[["subClonalRatio"]]))

#plot dNdS table
dNdStab <- read.csv(file="/Users/cross01/Desktop/dNdsTab.csv", header = TRUE, stringsAsFactors = FALSE)





