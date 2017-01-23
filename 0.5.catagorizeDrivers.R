# 1.inputs summary list of drivers and assigns driver mutation catagory

################### notes ###################
#SNVs and indels need to be merged and indels need to be manually inspected for duplicates

################# subroutines #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[1]])

#load driver gene list
annoDriverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table.05-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)

mutListDir <- "9.driverEvents/"
mutFile <- "lynch.somaticVariants.txt"

#assign mutation catagories to
mutName <- paste(sampleList[1,6], mutListDir, mutFile, sep="")
mutListTab <- read.table(file=mutName, header=TRUE, stringsAsFactors = FALSE, sep="\t")
mutListTab[10:12] <- NA
names(mutListTab)[10:12] <- c("COSMICfreq", "cat", "SNP-indel")

#SNP and common indel files, current verion: human_9606_b142_GRCh37p13
SNPfile <- "/Users/cross01/PhD/ReferenceGenome/00-All.txt"

for(currMut in 1:nrow(mutListTab)){
  trunkatingFlag <- 0
  currInfo <- mutListTab[currMut, "gene"]
  currInfo <- strsplit(currInfo, ",")
  if(length(currInfo[[1]]) > 1){
    currProt <- c()
    for(currAdd in 1:length(currInfo[[1]])){
      currProt[currAdd] <- strsplit(currInfo[[1]][currAdd], split = ":")[[1]][5]
    }
    if(length(unique(currProt)) == 1){
      currProt <- currProt[1]
    }
  }else{
    currProt <- strsplit(currInfo[[1]][1], split = ":")[[1]][5]
  }
  currGene <- strsplit(mutListTab[currMut, "gene"], ":")[[1]][1]
  currChom <- mutListTab[currMut, "chrom"]
  currPos <- mutListTab[currMut, "pos"]
  currAlt <- mutListTab[currMut, "alt"]
  
  #if X character present change to *
  for(currAss in 1:length(currProt)){
    if(gregexpr("X", currProt[currAss])[[1]][1] != -1){
      currProt <- gsub("X", "*", currProt, fixed=TRUE)
      trunkatingFlag <- 1
    }
  }
  
  #get mutations list for this gene
  mutListsName <- paste(sampleList[1,6], mutListDir, "mutationLists/", currGene, ".csv", sep="")
  if(file.exists(mutListsName)){
    MutList <- read.csv(file=mutListsName, header=TRUE, stringsAsFactors = FALSE)
  }else{
    mutListTab[currMut, "cat"] <- 3
    mutListTab[currMut, "COSMICfreq"] <- 0
    next
  }
  #remove * from mutation lists
  for(currRem in 1:nrow(MutList)){
    MutList[currRem, "AA.Mutation"] <- strsplit(MutList[currRem, "AA.Mutation"], split = "\\*")[[1]][1]
  }
  
  #get mutational frequencies
  totalMuts <- sum(MutList[["Count"]])
  MutList["Count"] <- MutList[["Count"]] / totalMuts
  
  #get current AA positions
  currAAlist <- c()
  for(currCh in 1:length(currProt)){
    AAflag <- 0
    currAA <- as.character()
    for(currChar in 1:nchar(currProt[currCh])){
      if(suppressWarnings(all(!is.na(as.numeric(substr(currProt[currCh], currChar, currChar)))))){
        currAA <- paste(currAA, substr(currProt[currCh], currChar, currChar), sep="")
        AAflag <- 1
      }else if(AAflag == 1){
        break
      }
    }
    currAAlist[currCh] <- currAA
  }
  
  
  #check if mutation is in mut list
  for(currCheck in 1:length(currAAlist)){
    AAflag <- 0
    if(currProt[currCheck] %in% MutList[["AA.Mutation"]]){
      mutListTab[currMut, "COSMICfreq"] <- max(MutList[MutList[["AA.Mutation"]] == currProt[currCheck], "Count"])
      mutListTab[currMut, "cat"] <- 1
      break
    }else if(trunkatingFlag == 1 | currAAlist[currCheck] %in% MutList[["Position"]]){
      mutListTab[currMut, "cat"] <- 2
      mutListTab[currMut, "COSMICfreq"] <- 0
      AAflag <- 1
    }else if(AAflag == 0){
      mutListTab[currMut, "cat"] <- 3
      mutListTab[currMut, "COSMICfreq"] <- 0
    }
  }
  
  
  #check if variant is a SNV or common indel
  print(paste("#### getting variants ", currChom, ":", currPos, " for set ", mutListTab[currMut, "Set"], " ####", sep=""))
  testStr <- paste('grep -w "', currChom, "\t", currPos, '" ', SNPfile, sep="")
  matchPos <- suppressWarnings(system(command = testStr, intern = TRUE))
  if(is.null(attributes(matchPos))){
    if(length(matchPos) > 1){
      matchTab <- data.frame(matrix(NA, ncol=6, nrow=length(matchPos)))
      names(matchTab) <- c("chrom", "pos", "rs", "ref", "alt", "stat")
      for(currMat in 1:nrow(matchTab)){
        matchTab[currMat, ] <- strsplit(matchPos[currMat], split = "\t")[[1]][1:6]
      }
    }else{
      matchTab <- data.frame(matrix(NA, ncol=6, nrow=1))
      names(matchTab) <- c("chrom", "pos", "rs", "ref", "alt", "stat")
      matchTab[1, ] <- strsplit(matchPos[1], split = "\t")[[1]][1:6]
    }
    
    #is this a SNP?
    if(mutListTab[currMut,1] == "nonsynonymous SNV" | mutListTab[currMut,1] == "stopgain" | mutListTab[currMut,1] == "synonymous SNV"){
      #it is an SNV
      currType <- "SNV"
    }else{
      currType <- "indel"
    }
    
    #determine if variant matches
    for(matchRow in 1:nrow(matchTab)){
      if(currType == "SNV" & nchar(matchTab[matchRow, "ref"]) == 1 & nchar(matchTab[matchRow, "alt"]) == 1){
        if(as.numeric(matchTab[matchRow, "chrom"]) == as.numeric(currChom) & as.numeric(matchTab[matchRow, "pos"]) == as.numeric(currPos) & matchTab[matchRow, "alt"] == mutListTab[currMut, "alt"]){
          mutListTab[currMut, "SNP-indel"] <- matchTab[matchRow, "rs"]
          break
        }
      }else if(currType == "indel" & matchTab[matchRow, "ref"] == mutListTab[currMut, "ref"] & (nchar(matchTab[matchRow, "alt"]) > 1 | matchTab[matchRow, "alt"] == "-")){
        mutListTab[currMut, "SNP-indel"] <- matchTab[matchRow, "rs"]
      }else{
        mutListTab[currMut, "SNP-indel"] <- paste(matchTab[matchRow, "rs"], "(SNP)", sep="")
      }
    }
  }else{
    mutListTab[currMut, "SNP-indel"] <- 0
  }
  
}  
 
#write output
mutFileTotal <- paste(sampleList[1,6], mutListDir, "lynch.somaticVariants.anno.txt", sep="")
write.table(mutListTab, file=mutFileTotal, sep="\t", quote = FALSE, row.names=FALSE, col.names = TRUE)

  