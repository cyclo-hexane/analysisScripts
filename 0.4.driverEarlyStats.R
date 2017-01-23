# uses txt vcf to determine whether driver genes are enriched on trunks by chance

######################## notes ########################
#
#
#
#
######################## libraries ########################


######################## subroutines ########################


######################## main program ########################


#input sampleList from commandline arguments
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
holdingDir <- "1.platypusCalls/exome.0.01/"
namePrepended <- ".snv.annoVar.exonic_variant_function.0.01.txt"
plotDir <- "1.platypusCalls/exome.0.01/"

sampleNames <- unique(sampleList[[1]])
#sampleNames <- sampleNames[-c(17:23)]
cancerNames <- sampleNames[c(1:3,5:11)]
lynchNames <- sampleNames[c(4,21:24)]
adenomaNames <- sampleNames[c(12:20)] 

######################## assess driver stats ########################


#read in driver gene list
fileName <- paste("Dropbox/MSeq-CRC/mseq.writeup/suppTables.final/table.S4-driverMutations.v2.csv", sep="")
driverMutList <- read.csv(file=fileName, header=TRUE, stringsAsFactors = FALSE)
driverMutList <- driverMutList[driverMutList[["cat"]]!="LOW", ]

#get gene names
for(currRow in 1:nrow(driverMutList)){
  driverMutList[currRow, 3] <- strsplit(driverMutList[currRow, 3], split = ":")[[1]][1] 
}


#excluding 2 biopsie sets
adenomaSet <- driverMutList[c(1:43), ]
carinomaSet <- driverMutList[c(44:51,59:92), ]
lynchSet <- driverMutList[c(52:58, 93:148), ]

#high confidence carcinomas
#carinomaSet <- carinomaSet[ (carinomaSet[["type"]]=="stopgain" | carinomaSet[["type"]]=="nonsynonymous SNV"), ]
#carinomaSet <- carinomaSet[carinomaSet[["COSMICfreq"]]!="A", ]
#numberCaDrivers <- nrow(carinomaSet[ carinomaSet[["cat"]]=="HIGH", ])
numberCaDrivers <- nrow(carinomaSet)
#noClonalCaDrivers <- nrow(carinomaSet[ carinomaSet[["phyloLoc"]]=="T" & carinomaSet[["cat"]]=="HIGH", ])
noClonalCaDrivers <- nrow(carinomaSet[ carinomaSet[["phyloLoc"]]=="T", ])

#high confidence lynch
#lynchSet <- lynchSet[ (lynchSet[["type"]]=="stopgain" | lynchSet[["type"]]=="nonsynonymous SNV"), ]
#lynchSet <- lynchSet[lynchSet[["COSMICfreq"]]!="A", ]
numberLyDrivers <- nrow(lynchSet[ lynchSet[["phyloLoc"]]!="GERMLINE", ])
#numberLyDrivers <- nrow(lynchSet)
noClonalLyDrivers <- nrow(lynchSet[ lynchSet[["phyloLoc"]]=="T", ])
#noClonalLyDrivers <- nrow(lynchSet[ lynchSet[["phyloLoc"]]=="T", ])

#high confidence
#adenomaSet <- adenomaSet[ (adenomaSet[["type"]]=="stopgain" | adenomaSet[["type"]]=="nonsynonymous SNV"), ]
#adenomaSet <- adenomaSet[adenomaSet[["COSMICfreq"]]!="A", ]
#numberAdDrivers <- nrow(adenomaSet[ adenomaSet[["cat"]]=="HIGH", ])
numberAdDrivers <- nrow(adenomaSet)
#noClonalAdDrivers <- nrow(adenomaSet[ adenomaSet[["phyloLoc"]]=="T" & adenomaSet[["cat"]]=="HIGH", ])
noClonalAdDrivers <- nrow(adenomaSet[ adenomaSet[["phyloLoc"]]=="T", ])



######## carcinoma analysis ###########

noCaTrunk <- 0
noCaBranchLeaf <- 0

#calculate driver stats for cancers
for(j in 1:length(cancerNames)){
	print(paste("#### analyzing sample ", cancerNames[j], " ####",sep=""))
	
	#subset main list
	subSample <- subset(sampleList, sampleList[1]==cancerNames[j])
	
	setName <- unique(subSample[[1]])
	noSamples <- subSample[1,8]
	normalIndex <- 1+subSample[1,7]
	samNames <- subSample[[2]]
	normalName <- samNames[normalIndex]
	samNamesNoNorm <- samNames[-normalIndex]
	
	#setup input/output names
	dataIn <- read.table(file=paste(subSample[1,6], holdingDir, setName,"/", setName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
	
	#name columns
  names(dataIn) <- c("line", "type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(samNames, ".NR", sep=""), samNames)
  dataIn[ncol(dataIn)+1] <- NA
  names(dataIn)[ncol(dataIn)] <- "phyloLoc"
  
  #assess phylogenetic locations and remove blank rows
  for(currRow in 1:nrow(dataIn)){
    assessRow <- dataIn[currRow, samNamesNoNorm]
    assessRow <- table(assessRow > 0) 
    if("TRUE" %in% names(assessRow)){
      if(as.numeric(assessRow["TRUE"]) == length(samNamesNoNorm)){
        #truncal
        dataIn[currRow, "phyloLoc"] <- "T"
      }else if(as.numeric(assessRow["TRUE"]) == 1){
        #leaf
        dataIn[currRow, "phyloLoc"] <- "L"
      }else{
        #branch
        dataIn[currRow, "phyloLoc"] <- "B"
      }
    }else{
    dataIn[currRow, "phyloLoc"] <- "R"
    } 
  }
  dataIn <- dataIn[dataIn[["phyloLoc"]] != "R", ]
  
  #get number of variants
  noCaTrunk <- noCaTrunk + nrow(dataIn[dataIn[["phyloLoc"]] == "T", ])
  noCaBranchLeaf <- noCaBranchLeaf + nrow(dataIn[dataIn[["phyloLoc"]] == "B" | dataIn[["phyloLoc"]] == "L", ])
}

#calculate stats
binomPvalue <- binom.test(x=noClonalCaDrivers, n=numberCaDrivers, p =(noCaTrunk/(noCaTrunk+noCaBranchLeaf)), alternative = "two.sided")






######## Lynch syndrome analysis ##########

noLyTrunk <- 0
noLyBranchLeaf <- 0

#calculate driver stats for cancers
for(j in 1:length(lynchNames)){
  print(paste("#### analyzing sample ", lynchNames[j], " ####",sep=""))
  
  #subset main list
  subSample <- subset(sampleList, sampleList[1]==lynchNames[j])
  
  setName <- unique(subSample[[1]])
  noSamples <- subSample[1,8]
  normalIndex <- 1+subSample[1,7]
  samNames <- subSample[[2]]
  normalName <- samNames[normalIndex]
  samNamesNoNorm <- samNames[-normalIndex]
  
  #setup input/output names
  dataIn <- read.table(file=paste(subSample[1,6], holdingDir, setName,"/", setName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
  
  #name columns
  names(dataIn) <- c("line", "type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(samNames, ".NR", sep=""), samNames)
  dataIn[ncol(dataIn)+1] <- NA
  names(dataIn)[ncol(dataIn)] <- "phyloLoc"
  
  #assess phylogenetic locations and remove blank rows
  for(currRow in 1:nrow(dataIn)){
    assessRow <- dataIn[currRow, samNamesNoNorm]
    assessRow <- table(assessRow > 0) 
    if("TRUE" %in% names(assessRow)){
      if(as.numeric(assessRow["TRUE"]) == length(samNamesNoNorm)){
        #truncal
        dataIn[currRow, "phyloLoc"] <- "T"
      }else if(as.numeric(assessRow["TRUE"]) == 1){
        #leaf
        dataIn[currRow, "phyloLoc"] <- "L"
      }else{
        #branch
        dataIn[currRow, "phyloLoc"] <- "B"
      }
    }else{
      dataIn[currRow, "phyloLoc"] <- "R"
    } 
  }
  dataIn <- dataIn[dataIn[["phyloLoc"]] != "R", ]
  
  #get number of variants
  noLyTrunk <- noCaTrunk + nrow(dataIn[dataIn[["phyloLoc"]] == "T", ])
  noLyBranchLeaf <- noCaBranchLeaf + nrow(dataIn[dataIn[["phyloLoc"]] == "B" | dataIn[["phyloLoc"]] == "L", ])
}

#calculate stats
binomPvalue <- binom.test(x=noClonalLyDrivers, n=numberLyDrivers, p =(noLyTrunk/(noLyTrunk+noLyBranchLeaf)), alternative = "two.sided")






######### adenoma analysis ###############

noAdTrunk <- 0
noAdBranchLeaf <- 0

for(j in 1:length(adenomaNames)){
  print(paste("#### analyzing sample ", adenomaNames[j], " ####",sep=""))
  
  #subset main list
  subSample <- subset(sampleList, sampleList[1]==adenomaNames[j])
  
  setName <- unique(subSample[[1]])
  noSamples <- subSample[1,8]
  normalIndex <- 1+subSample[1,7]
  samNames <- subSample[[2]]
  normalName <- samNames[normalIndex]
  samNamesNoNorm <- samNames[-normalIndex]
  
  #setup input/output names
  dataIn <- read.table(file=paste(subSample[1,6], holdingDir, setName,"/", setName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
  
  #name columns
  names(dataIn) <- c("line", "type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(samNames, ".NR", sep=""), samNames)
  dataIn[ncol(dataIn)+1] <- NA
  names(dataIn)[ncol(dataIn)] <- "phyloLoc"
  
  #assess phylogenetic locations and remove blank rows
  for(currRow in 1:nrow(dataIn)){
    assessRow <- dataIn[currRow, samNamesNoNorm]
    assessRow <- table(assessRow > 0) 
    if("TRUE" %in% names(assessRow)){
      if(as.numeric(assessRow["TRUE"]) == length(samNamesNoNorm)){
        #truncal
        dataIn[currRow, "phyloLoc"] <- "T"
      }else if(as.numeric(assessRow["TRUE"]) == 1){
        #leaf
        dataIn[currRow, "phyloLoc"] <- "B"
      }else{
        #branch
        dataIn[currRow, "phyloLoc"] <- "L"
      }
    }else{
      dataIn[currRow, "phyloLoc"] <- "R"
    } 
  }
  dataIn <- dataIn[dataIn[["phyloLoc"]] != "R", ]
  
  #get number of variants
  noAdTrunk <- noAdTrunk + nrow(dataIn[dataIn[["phyloLoc"]] == "T", ])
  noAdBranchLeaf <- noAdBranchLeaf + nrow(dataIn[dataIn[["phyloLoc"]] == "B" | dataIn[["phyloLoc"]] == "L", ])
}

#calculate stats
binomPvalue <- binom.test(x=noClonalAdDrivers, n=numberAdDrivers, p=(noAdTrunk/(noAdTrunk+noAdBranchLeaf)), alternative = "two.sided")

