# takes list of driver mutations and subsets vcf file
# outputs counts of variants and clonal proportions
# produces a heatmap table of canonical mutations

################### notes ###################
#
# uses an indel, SNV and CNV table and driver list to produce summary table
#
#
#
################# subroutines #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[1]])

#load driver gene list
annoDriverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table5-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)

SNVdir <- "9.driverEvents/SNVs/"
SNVPrep <- ".snv.Platypus.drivers.txt"

indelDir <- "9.driverEvents/indels/"
indelPrep <- ".combIndels.txt"

CNVdir <- "8.CNVphylogenetics/CNVphylo/"
CNVprep <- ".phyloCNVs.anno.csv"

countsTable <- data.frame(matrix(NA, ncol=7, nrow=length(sampleNames)))
names(countsTable) <- c("total", "truncal", "branch", "leaf", "noDrivers", "APC", "clonalProp")
row.names(countsTable) <- sampleNames

for(currSet in 1:length(sampleNames)){
  currentSample <- subset(sampleList, sampleList[1]==sampleNames[currSet])
  nonNormalNames <- currentSample[[2]]
  nonNormalNames <- nonNormalNames[nonNormalNames!=nonNormalNames[(currentSample[1,7]+1)]]
  
  #get SNVs
  SNVfileName <- paste(currentSample[1,6], SNVdir, currentSample[1,1], SNVPrep, sep="")
  SNVfile <- read.table(file=SNVfileName, header=FALSE, stringsAsFactors = FALSE, sep="\t", colClasses = "character")
  names(SNVfile) <- c("varType", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(currentSample[[2]], ".NR", sep=""), currentSample[[2]])
  
  #total variants
  countsTable[currSet, 1] <- nrow(SNVfile)
  
  #get gene names
  keepRow <- c()
  keepCounter <- 1
  APCflag <- 0
  for(currRow in 1:nrow(SNVfile)){
    SNVfile[currRow, 3] <- strsplit(SNVfile[currRow, 3], split=":")[[1]][1]
    if(SNVfile[currRow, 3] %in% annoDriverList[[1]]){
      keepRow[keepCounter] <- currRow
      keepCounter <- keepCounter + 1
      
      #mark if APC
      if(SNVfile[currRow, 3] == "APC"){
        APCflag <- 1
      }
    }
  }
  driverFile <- SNVfile[keepRow, ]
  driverFile <- driverFile[-1]
  driverFileName <- paste(currentSample[1,6], SNVdir, currentSample[1,1], "/", currentSample[1,1], ".snv.annoVar.exonic_variant_function.drivers.txt", sep="")
  #write.table(driverFile, file=driverFileName, sep="\t", quote = FALSE, row.names=FALSE, col.names = FALSE)
  
  #number of drivers
  countsTable[currSet, 5] <- nrow(driverFile)
  
  #APC present?
  if(APCflag == 1){
    countsTable[currSet, 6] <- "yes"
  }else{
    countsTable[currSet, 6] <- "no"
  }
  
  #get phylogenetic locations
  SNVfile[ncol(SNVfile)+1] <- NA
  names(SNVfile)[ncol(SNVfile)] <- "phyloLoc"
  for(currVar in 1:nrow(SNVfile)){
    assessRow <- table(SNVfile[currVar, nonNormalNames] > 0)
    if("FALSE" %in% names(assessRow)){
      if(assessRow["TRUE"]==1){
        SNVfile[currVar, "phyloLoc"] <- "L"
      }else{
        SNVfile[currVar, "phyloLoc"] <- "B"
      }
    }else{
      SNVfile[currVar, "phyloLoc"] <- "T"
    }
  }
  
  #number of trunk branch and leaves
  countsTable[currSet, 2] <- nrow(SNVfile[SNVfile[["phyloLoc"]]=="T", ])
  countsTable[currSet, 3] <- nrow(SNVfile[SNVfile[["phyloLoc"]]=="B", ])
  countsTable[currSet, 4] <- nrow(SNVfile[SNVfile[["phyloLoc"]]=="L", ])
  
  #get proportion
  countsTable[currSet, 7] <- countsTable[currSet, 2] / countsTable[currSet, 1]
}

#write table
fileOut <- paste(currentSample[1,6], "2.phylogenetics/clonalProportions.total.csv", sep="")
write.csv(countsTable, file=fileOut, row.names=TRUE)


#make pathway table for each sample where
mutTabList <- data.frame(matrix(NA, nrow=0, ncol=nrow(annoDriverList)))

for(currSam in 1:length(sampleNames)){
  currentSample <- subset(sampleList, sampleList[1]==sampleNames[currSam])
  
  varTable <- data.frame(matrix(NA, ncol=nrow(annoDriverList), nrow=5))
  colnames(varTable) <- annoDriverList[[1]]
  rownames(varTable) <- paste(sampleNames[currSam], c("SNV1", "SNV2", "indel1", "indel2", "CNV"), sep=".")
  
  #get variant calls
  SNVfile <- paste(currentSample[1,6], SNVdir, currentSample[1,1], SNVPrep, sep="")
  SNVdata <- read.table(file=SNVfile, header=FALSE, stringsAsFactors = FALSE, sep="\t")
  
  SNVdriverList <- c()
  #get SNV driver names
  for(currDri in 1:nrow(SNVdata)){
    SNVdriverList[currDri] <- strsplit(SNVdata[currDri, 2], split = ":")[[1]][1]
  }
  
  #get indels
  indelFile <- paste(currentSample[1,6], indelDir, currentSample[1,1], indelPrep, sep="")
  if(file.exists(indelFile)){
    indelData <- read.table(file=indelFile, header=FALSE, stringsAsFactors = FALSE, sep="\t")
  }else{
    indelData <- NA
  }

  indelDriverList <- c()
  #get indel driver names
  if(!is.null(nrow(indelData))){
    for(currDri in 1:nrow(indelData)){
      indelDriverList[currDri] <- strsplit(indelData[currDri, 2], split = ":")[[1]][1]
    }
  }
  
  CNVfile <- paste(currentSample[1,6], CNVdir, currentSample[1,1], CNVprep, sep="")
  CNVdata <- read.table(file=CNVfile, header=FALSE, stringsAsFactors = FALSE, sep=",")
  CNVdata <- CNVdata[CNVdata[(ncol(CNVdata)-1)]!="none", ]
  
  #populate table with variants
  for(currMut in 1:ncol(varTable)){
    tempMut <- names(varTable)[currMut]
    
    #check if mutation in SNV data
    if(tempMut %in% SNVdriverList){
      rowIndex <- which(tempMut == SNVdriverList)
      tempSNV <- SNVdata[rowIndex, ]
      
      #get HGV nomenclature
      for(currRow in 1:nrow(tempSNV)){
        tempSNV[currRow, 3] <- paste(tempSNV[currRow, 1], ":", strsplit(tempSNV[currRow, 2], split = ":")[[1]][4], sep="")
      }
      
      if(nrow(tempSNV) == 1){
        varTable[1, currMut] <- tempSNV[1,3] 
      }else{
        varTable[c(1:2), currMut] <- tempSNV[c(1:2),3]
      } 
    }
    
    #check if mutation in indel list
    if(!is.null(nrow(indelData))){
      if(tempMut %in% indelDriverList){
        rowIndex <- which(tempMut == indelDriverList)
        tempIndel <- indelData[rowIndex, ]
        
        #get HGV nomenclature
        for(currRow in 1:nrow(tempIndel)){
          tempIndel[currRow, 3] <- paste(tempIndel[currRow, 1], ":", strsplit(tempIndel[currRow, 2], split = ":")[[1]][3], sep="")
        }
        
        if(nrow(tempIndel) == 1){
          varTable[3, currMut] <- tempIndel[1, 3]
        }else{
          varTable[c(3:4), currMut] <- tempIndel[c(1:2), 3]
        } 
      }
    }
    
    #check for CNV status
    for(currCNV in 1:nrow(CNVdata)){
      #parse gene names
      geneNameList <- strsplit(CNVdata[currCNV, (ncol(CNVdata)-1)], split = ":")[[1]]
      
      if(tempMut %in% geneNameList){
        if(CNVdata[currCNV, 6] == 0 & CNVdata[currCNV, 7] == 1){
          CNVtype <- "LOH"
        }else if(CNVdata[currCNV, 6] == 0 & CNVdata[currCNV, 7] > 1){
          CNVtype <- "BLOH"
        }else{
          CNVtype <- "gain"
        }
        varTable[5, currMut] <- CNVtype
        break
      }
    }
  }
  
  #save table to list
  mutTabList <- rbind(mutTabList, varTable)

  #remove mutations not represented
  remList <- c()
  remCounter <- 1
  for(currMut in 1:ncol(mutTabList)){
    assessRow <- table(is.na(mutTabList[[currMut]]))
    if(FALSE %in% names(assessRow)){
      #do nothing
    }else{
      remList[remCounter] <- currMut
      remCounter <- remCounter + 1
    }
  }
  
  #write table
  fileOut <- paste(currentSample[1,6], "9.driverEvents/driverRepresentationList.txt", sep="")
  write.table(mutTabList, file=fileOut, sep="\t", quote = FALSE, col.names = TRUE, row.names=TRUE)
}



