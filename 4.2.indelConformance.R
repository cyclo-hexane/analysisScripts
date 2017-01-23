# produces a conformance list of variants and stats

################### notes ###################
#
# call variants using Platypus and Scalpel (perl scripts)
#
# 4.0.indelPlatypusAnnotation (correct columns and annotate using annovar)
#       |
#       v
# 4.1.indelCompileScalpel (compileTables for Scalpel calls)
#       |
#       v
# 4.2.indelConformance (compare Platypus to Scalpel and get conformance table, combine driver gene lists) 
#       |
#       v
# 4.3.phylogeneticsPrep.indels (make trees using PAUP)
#
#
################# libraries #################



################# subroutines #################



############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.indels.csv", header=FALSE, stringsAsFactors=FALSE)

sampleNames <- unique(sampleList[[1]])
#sampleNames <- sampleNames[c(2:5, 13:16)]

#indel platypus names
platDir <- "1.platypusCalls/indels/"
platName <- ".exonic_variant_function.txt"
platDrivers <- ".indel.Platypus.drivers.txt"

#indel scalpel names
scalpelDir <- "1d.scalpelCalls/"
scalpelName <- ".indel.Scalpel.exome.txt"
scalpelDrivers <- ".indel.Scalpel.drivers.txt"

conformanceDir <- "9.driverEvents/indels/"

conformanceData <- data.frame(matrix(NA, nrow=length(sampleNames), ncol=5))
names(conformanceData) <- c("platypusIndels", "scalpelIndels", "conformance", "propScalpel-Plat", "propPlat-Scalpel")
row.names(conformanceData) <- sampleNames

#compare call sets
for(currSam in 1:length(sampleNames)){
  print(paste("#### comparing indels for sample ", sampleNames[currSam], " ####",sep=""))
  
  subSample <- subset(sampleList, sampleList[1]==sampleNames[currSam])
  noSamples <- subSample[1,8]
  normalName <- subSample[(subSample[1,7]+1), 2]
  setNames <- subSample[[2]]
  setNames <- setNames[setNames!=normalName]
  
  #get current platypus indel vcf (.txt) file
  platIndelFile <- paste(subSample[1,6], platDir, sampleNames[currSam], "/", sampleNames[currSam], platName, sep="")
  platIndelData <- read.table(file=platIndelFile, header=TRUE, stringsAsFactors=FALSE, sep="\t")

  #get scalpel calls
  scalIndelFile <- paste(subSample[1,6], scalpelDir, sampleNames[currSam], "/", sampleNames[currSam], scalpelName, sep="")
  scalIndelData <- read.table(file=scalIndelFile, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  
  #compare files
  checkTable <- scalIndelData[c(1:8, ncol(scalIndelData))]
  checkTable[ncol(checkTable)+1] <- NA
  names(checkTable)[ncol(checkTable)] <- "present"
  
  for(currCheck in 1:nrow(checkTable)){
    currChrom <- checkTable[currCheck, "chrom"]
    currStart <- checkTable[currCheck, "posStart"]
    currEnd <- checkTable[currCheck, "posEnd"]
    currGene <- checkTable[currCheck, "gene"]
    currGene <- strsplit(currGene, ":")[[1]][1]
    currType <- checkTable[currCheck, "indel"]
    
    #look for current indel in platypus calls
    platSub <- subset(platIndelData, platIndelData["chrom"]==currChrom)
    if(currStart %in% platSub[["posStart"]] | currStart %in% (platSub[["posStart"]]-1) | currStart %in% (platSub[["posStart"]]+1) | currStart %in% (platSub[["posStart"]]-2) | currStart %in% (platSub[["posStart"]]+2)){
      #check sample presence matches
      currentIndel <- scalIndelData[currCheck, paste(setNames, ".NV", sep="")]
      scalPresence <- setNames[currentIndel>0]
      
      #presence in platpus calls
      platPresence <- platSub[platSub[["posStart"]]==currStart | (platSub[["posStart"]]-1)==currStart | (platSub[["posStart"]]-2)==currStart | (platSub[["posStart"]]+1)==currStart | (platSub[["posStart"]]+2)==currStart, paste(setNames, ".NV", sep="")]
      platPresence <- setNames[platPresence>0]

      #if sample presence match
      if(identical(platPresence, scalPresence)){
        checkTable[currCheck, "present"] <- 2
      }else{
        #mark as partial match
        checkTable[currCheck, "present"] <- 1
      }
      
    }else{
      checkTable[currCheck, "present"] <- 0
    } 
  }
  
  #get stats
  conformanceData[currSam, 1] <- nrow(platIndelData)
  conformanceData[currSam, 2] <- nrow(scalIndelData)
  conformanceData[currSam, 3] <- nrow(checkTable[checkTable[["present"]]>0,])
  conformanceData[currSam, 4] <- conformanceData[currSam, 3] / nrow(scalIndelData)
  conformanceData[currSam, 5] <- conformanceData[currSam, 3] / nrow(platIndelData)
  
  #get conformance Data
  #confScalpelData <- scalIndelData[checkTable[["present"]]>0,]
}

#output new table
comfIndelFile <- paste(subSample[1,6], conformanceDir, "indel.conformance.txt", sep="")
write.table(conformanceData, file=comfIndelFile, sep="\t", quote = FALSE, col.names=TRUE, row.names=TRUE)


#combine driver mutation files
for(combDriver in 1:length(sampleNames)){
  subSample <- subset(sampleList, sampleList[1]==sampleNames[combDriver])
  
  combinedFile <- paste(subSample[1,6], conformanceDir, subSample[1,1], ".combIndels.txt", sep="")
  platFile <- paste(subSample[1,6], platDir, subSample[1,1], "/", subSample[1,1], platDrivers, sep="")
  scalpelFile <- paste(subSample[1,6], scalpelDir, subSample[1,1], "/",  subSample[1,1], scalpelDrivers, sep="")
  
  if(!file.exists(platFile)){
    if(!file.exists(scalpelFile)){
      next()
    }else{
      scalpelIn <- read.table(file=scalpelFile, sep="\t", header=TRUE, stringsAsFactors = FALSE, colClasses = c(rep("character", 3), rep("numeric", 2), rep("character", 3), rep("numeric", (length(setNames)*2))) )
      write.table(scalpelIn, file=combinedFile, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
      next()
    }
  }else if(!file.exists(scalpelFile)){
    platIn <- read.table(file=platFile, sep="\t", header=TRUE, stringsAsFactors = FALSE,colClasses = c(rep("character", 3), rep("numeric", 2), rep("character", 3), rep("numeric", (length(setNames)*2))))
    write.table(platIn, file=combinedFile, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    next()
  }else{
    scalpelIn <- read.table(file=scalpelFile, sep="\t", header=TRUE, stringsAsFactors = FALSE, colClasses = c(rep("character", 3), rep("numeric", 2), rep("character", 3), rep("numeric", (length(setNames)*2))))
    platIn <- read.table(file=platFile, sep="\t", header=TRUE, stringsAsFactors = FALSE, colClasses = c(rep("character", 3), rep("numeric", 2), rep("character", 3), rep("numeric", (length(setNames)*2))))
    
    #merge tables
    combinedTab <- rbind(scalpelIn, platIn)
    combinedTab <- combinedTab[order(combinedTab[["chrom"]], combinedTab[["posStart"]]), ]
    write.table(combinedTab, file=combinedFile, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }  
}


