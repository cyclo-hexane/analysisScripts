# 1.takes list of driver mutations and subsets vcf file
# 2.also counts non-synonymous variants and outputs to table
# 3.inputs summary list of drivers and assigns driver mutation catagory

################### notes ###################


################# subroutines #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterAdinCaList.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/sampleList.Set.09.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[1]])
#sampleNames <- sampleNames[-15]

#load driver gene list
annoDriverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table5-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)

#Filedir <- "1.platypusCalls/indels/"
#FilePrep <- ".exonic_variant_function.txt"
#outPrep <- ".indel.Platypus.drivers.txt"

Filedir <- "1.platypusCalls/exome/"
FilePrep <- ".snv.annoVar.exonic_variant_function.txt"
outPrep <- ".driverSNVs.txt"

mutListDir <- "9.driverEvents/"
mutFile <- "driversIndelCounts.txt"

#Filedir <- "1d.scalpelCalls/"
#FilePrep <- ".indel.Scalpel.exome.txt"

geneNameCol <- 3

totalDrivers <- data.frame(matrix(NA, nrow=0, ncol=9))
names(totalDrivers) <- c("type", "gene", "chrom", "posStart", "posEnd", "ref", "alt", "phyloLoc", "sampleSet")

#non-synonymous variants
totalNSvars <- data.frame(matrix(NA, nrow=length(sampleNames), ncol=2))
names(totalNSvars) <- c("set", "nonSyn")

for(currSet in 1:length(sampleNames)){
  currentSample <- subset(sampleList, sampleList[1]==sampleNames[currSet])
  nonNorSamples <- currentSample[[2]]
  nonNorSamples <- nonNorSamples[-(currentSample[1,7]+1)]
  noSamples <- length(nonNorSamples)
  
  #get CNV and timings
  fileName <- paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], FilePrep, sep="")
  dataFile <- read.table(file=fileName, header=FALSE, stringsAsFactors = FALSE, sep="\t")
  
  totalNSvars[currSet, 1] <- currentSample[1,1]
  totalNSvars[currSet, 2] <- nrow(dataFile)
  
  #get gene names
  keepRow <- c()
  keepCounter <- 1
  for(currRow in 1:nrow(dataFile)){
    geneTemp <- strsplit(dataFile[currRow, geneNameCol], split=":")
    #dataFile[currRow, 2] <- paste(geneTemp[[1]][1], ":", geneTemp[[1]][3], ":", geneTemp[[1]][4], sep="")
    #dataFile[currRow, geneNameCol] <- strsplit(dataFile[currRow, geneNameCol], split=":")[[1]][5]
    if(geneTemp[[1]][1] %in% annoDriverList[[1]] | geneTemp[[1]][1] %in% annoDriverList[[2]]){
      keepRow[keepCounter] <- currRow
      keepCounter <- keepCounter + 1
    }
  }
  driverFile <- dataFile[keepRow, ]
  driverFile <- driverFile[-1]
  
  if(nrow(driverFile)!=0){
    driverFileName <- paste(currentSample[1,6], Filedir, currentSample[1,1], "/", currentSample[1,1], outPrep, sep="")
    write.table(driverFile, file=driverFileName, sep="\t", quote = FALSE, row.names=FALSE, col.names =FALSE)
  
    names(driverFile) <- c("type", "gene", "chrom", "pos", "pos2", "ref", "alt",  paste(currentSample[[2]], ".NV", sep=""), currentSample[[2]])
    
    #get driver name and phyloLoc
    for(currRow in 1:nrow(driverFile)){
      driverFile[currRow, 2] <- strsplit(driverFile[currRow, 2], split=":")[[1]][1]
      assessRow <-  driverFile[currRow, nonNorSamples]
      if(as.numeric(table(assessRow > 0)["TRUE"]) == noSamples){
        driverFile[currRow, 8] <- "T"
      }else if(as.numeric(table(assessRow > 0)["TRUE"]) == 1){
        driverFile[currRow, 8] <- "L"
      }else{
        driverFile[currRow, 8] <- "B"
      }
    }
    
    #save driver gene to main list
    driverFile[9] <- currentSample[1,1]
    names(driverFile)[8:9] <- c("phyloLoc", "sampleSet")
    totalDrivers <- rbind(totalDrivers, driverFile[1:9])
  }
}

#write total table
names(totalDrivers) <- c("type" ,"gene", "chrom", "pos", "pos2", "ref", "alt", "phyloLoc", "Set")
driverFileTotal <- paste(currentSample[1,6], Filedir, "Ha.drivers.txt", sep="")
write.table(totalDrivers, file=driverFileTotal, sep="\t", quote = FALSE, row.names=FALSE, col.names = TRUE)

#write NS file
nsFileTotal <- paste(currentSample[1,6], Filedir, "Ha.nsVars.txt", sep="")
write.table(totalNSvars, file=nsFileTotal, sep="\t", quote = FALSE, row.names=FALSE, col.names = TRUE)


#now assign mutation catagories to

mutName <- paste(sampleList[1,6], mutListDir, mutFile, sep="")
mutListTab <- read.table(file=mutName, header=TRUE, stringsAsFactors = FALSE, sep="\t")
mutListTab[10:11] <- NA
names(mutListTab)[10:11] <- c("COSMICfreq", "cat")
  
for(currMut in 1:nrow(mutListTab)){
  trunkFlag <- 0
  currGene <- mutListTab[currMut, "gene"]
  currProt <- strsplit(mutListTab[currMut, "protein"], split = ",")[[1]][1]
  
  #if X character present change to *
  if(gregexpr("X", currProt)[[1]][1] != -1 | gregexpr("fs", currProt)[[1]][1] != -1){
    currProt <- gsub("X", "*", currProt, fixed=TRUE)
    currProt <- gsub("fs", "*", currProt, fixed=TRUE)
    trunkFlag <- 1
  }
  
  #get mutations list for this gene
  mutListsName <- paste(sampleList[1,6], mutListDir, "mutationLists/", currGene, ".csv", sep="")
  MutList <- read.csv(file=mutListsName, header=TRUE, stringsAsFactors = FALSE)
  
  #get mutational frequencies
  totalMuts <- sum(MutList[["Count"]])
  MutList["Count"] <- MutList[["Count"]] / totalMuts
  
  #get current AA position
  currAA <- as.character()
  for(currChar in 1:nchar(currProt)){
    if(suppressWarnings(all(!is.na(as.numeric(substr(currProt, currChar, currChar)))))){
      currAA <- paste(currAA, substr(currProt, currChar, currChar), sep="")
    }
  }
  currAA <- as.numeric(currAA)
  
  #check if mutation is in mut list
  if(currProt %in% MutList[["AA.Mutation"]]){
    mutListTab[currMut, "COSMICfreq"] <- max(MutList[MutList[["AA.Mutation"]] == currProt, "Count"])
    mutListTab[currMut, "cat"] <- 1
  }else if(trunkFlag == 1 | currAA %in% MutList[["Position"]]){
    mutListTab[currMut, "cat"] <- 2
    mutListTab[currMut, "COSMICfreq"] <- 0
  }else{
    mutListTab[currMut, "cat"] <- 3
    mutListTab[currMut, "COSMICfreq"] <- 0
  }
}  
 
#write output
mutFileTotal <- paste(sampleList[1,6], mutListDir, "Ha.driversCats.SNVs.txt", sep="")
write.table(mutListTab, file=mutFileTotal, sep="\t", quote = FALSE, row.names=FALSE, col.names = TRUE)

  
  