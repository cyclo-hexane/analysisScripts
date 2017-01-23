# compiles a joint list for a given indel set

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
# 4.2.indelConformance (compare Platypus to Scalpel and get conformance table) 
#       |
#       v
# 4.3.phylogeneticsPrep.indels (make trees using PAUP)
#
#
################# libraries #################



################# subroutines #################



############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/archive/sampleList.Set.09.filtered.csv", header=FALSE, stringsAsFactors=FALSE)

sampleNames <- unique(sampleList[[1]])
#sampleNames <- sampleNames[24]

#driver gene list
driverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table.05-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)

#indel names
holdingDir <- "1d.scalpelCalls/"

#indelFileNames <- ".indel.annovar.variant_function.txt"
indelFileNames <- ".indel.annovar.exonic_variant_function.txt"

genomeFlag <- "exome"

#compile 
for(currSam in 1:length(sampleNames)){
  print(paste("#### compiling indel table for ", sampleNames[currSam], " ####",sep=""))
  
  subSample <- subset(sampleList, sampleList[1]==sampleNames[currSam])
  noSamples <- nrow(subSample)-1
  setNames <- subSample[[2]]
  normalIndex <- subSample[1, 7]+1
  setNames <- setNames[-normalIndex]
  
  #get files into list
  indelList <- as.list(NA)
  locationList <- data.frame(matrix(NA, nrow=0, ncol=8))
  names(locationList) <- c("type", "gene", "chrom", "posStart", "posEnd", "ref", "alt", "indel")
  for(currFile in 1:length(setNames)){
    #get current indel vcf (.txt) file
    indelFile <- paste(subSample[1,6], holdingDir, sampleNames[currSam], "/", setNames[currFile], indelFileNames, sep="")
    indelDataIn <- read.table(file=indelFile, header=FALSE, stringsAsFactors=FALSE, sep="\t")
    
    if(genomeFlag=="exome"){
      indelDataIn <- indelDataIn[-1]
    }
    
    #filter for normal variants
    indelDataIn <- indelDataIn[indelDataIn[[18]]=="no", ]
    
    #add to location
    addLoc <- indelDataIn[1:7]
    locationList <- rbind(indelDataIn[c(1:7, 10)], locationList)
    
    #add to list
    indelList[[currFile]] <- indelDataIn
  }
  #get mutation names
  #for(currMut in 1:nrow(locationList)){
  #  mutTemp <- strsplit(locationList[currMut, 2], ":")
  #  locationList[currMut, 2] <- paste(mutTemp[[1]][1], ":", mutTemp[[1]][3], ":", mutTemp[[1]][4], sep="")
  #}
  
  #get unique locations
  locationList <- unique(locationList)
  locationList <- locationList[order(locationList[[3]], locationList[[4]]), ]
  
  #merge table
  locationList[9:(8+(noSamples*2))] <- NA
  names(locationList) <- c("type", "gene", "chrom", "posStart", "posEnd", "ref", "alt", "indel", paste(setNames, ".NR", sep=""), paste(setNames, ".NV", sep=""))
  
  #populate location list
  for(currAdd in 1:length(setNames)){
    currentIndels <- indelList[[currAdd]]
    
    #add to location table
    for(currRow in 1:nrow(locationList)){
      currChrom <- locationList[currRow, 3]
      currPosStart <- locationList[currRow, 4]
      currPosEnd <- locationList[currRow, 5]
      
      #subset current indel list
      subIndels <- currentIndels[currentIndels[[3]]==currChrom, ]
      if(nrow(subIndels) == 0){
        locationList[currRow, paste(setNames[currAdd], ".NV", sep="")] <- 0
        locationList[currRow, paste(setNames[currAdd], ".NR", sep="")] <- 1
      }else if(currPosStart %in% subIndels[[4]]){
        foundVar <- subIndels[subIndels[[4]]==currPosStart, ]
        
        #check it is correct variant
        if(locationList[currRow, 6]==foundVar[1,6] & locationList[currRow, 7]==foundVar[1,7]){
          #mark variant VAF details
          locationList[currRow, paste(setNames[currAdd], ".NV", sep="")] <- foundVar[1, 12]
          locationList[currRow, paste(setNames[currAdd], ".NR", sep="")] <- foundVar[1, 12] + foundVar[1, 14]
        }else{
          #mark as absent
          locationList[currRow, paste(setNames[currAdd], ".NV", sep="")] <- 0
          locationList[currRow, paste(setNames[currAdd], ".NR", sep="")] <- foundVar[1, 12] + foundVar[1, 14]
        }
      }else{
        #mark as absent
        locationList[currRow, paste(setNames[currAdd], ".NV", sep="")] <- 0
        locationList[currRow, paste(setNames[currAdd], ".NR", sep="")] <- 1
      }
    }
  }
  
  #remove rows with zero total, these are indels reported at the same location within one sample 
  #remList <- c()
  #remCounter <- 1
  #for(currRem in 1:nrow(locationList)){
  #  assessRow <- sum(locationList[currRem, setNames])
  #  if(assessRow==0){
  #    remList[remCounter] <- currRem
  #    remCounter <- remCounter + 1
  #  }
  #}
  #locationList <- locationList[-remList, ]
  
  #move indel column to last column
  #locationList[ncol(locationList)+1] <- locationList[["indel"]]
  #names(locationList)[ncol(locationList)] <- "indelType"
  #remCol <- which("indel" == names(locationList))
  #locationList <- locationList[-remCol]
  
  #output new table
  #indelOutFile <- paste(subSample[1,6], holdingDir, sampleNames[currSam], "/", sampleNames[currSam], ".indel.Scalpel.txt", sep="")
  indelOutFile <- paste(subSample[1,6], holdingDir, sampleNames[currSam], "/", sampleNames[currSam], ".indel.Scalpel.exome.txt", sep="")
  write.table(locationList, file=indelOutFile, sep="\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
  
  #subset variants for driver genes only
  keepList <- c()
  keepCounter <- 1
  for(currIndel in 1:nrow(locationList)){
    currentGene <- strsplit(locationList[currIndel, "gene"], ":")[[1]][1]
    if(currentGene %in% driverList[[1]]){
      keepList[keepCounter] <- currIndel
      keepCounter <- keepCounter + 1
    }
  }
  locationList <- locationList[keepList, ]
  
  #get unique drivers
  locationList <- locationList[rownames(locationList) %in% rownames(unique(locationList[3:8])), ]
  
  if(nrow(locationList) != 0){
    indelDriverFile <- paste(subSample[1,6], holdingDir, sampleNames[currSam], "/", sampleNames[currSam], ".indel.Scalpel.drivers.txt", sep="")
    write.table(locationList, file=indelDriverFile, sep="\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
    
    #now format table to enable catagorization
    locationList <- locationList[1:9]
    names(locationList) <- c("gene", "change", "protein", "chrom", "pos", "exon", "ref", "alt", "set")
    locationList[9] <- sampleNames[currSam]
    locationList[8] <- locationList[7]
    locationList[7] <- locationList[6]
    locationList[4] <- locationList[3]
    locationList[3] <- locationList[2]
    locationList[2] <- locationList[1]
    for(currMut in 1:nrow(locationList)){
      locationList[currMut, "gene"] <- strsplit(locationList[currMut, "protein"], split = ":")[[1]][1]
      locationList[currMut, "exon"] <- strsplit(locationList[currMut, "protein"], split = ":")[[1]][3]
      locationList[currMut, "protein"] <- strsplit(locationList[currMut, "protein"], split = ":")[[1]][5]
    }
    locationList <- unique(locationList)
    
    indelDriverFile <- paste(subSample[1,6], holdingDir, sampleNames[currSam], "/", sampleNames[currSam], ".driversIndelCounts.txt", sep="")
    write.table(locationList, file=indelDriverFile, sep="\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
  }
}


