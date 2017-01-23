#script takes lynch polyp sample list of vcfs and removes variants also found in cancer

######################## notes ########################
#
#
######################## libraries ########################


######################## subroutines ########################


######################## main program ########################


sampleList <- read.csv(file="~/PhD/CRCproject/masterLynchList.polyps.filt.csv", header=FALSE, stringsAsFactors=FALSE)

setNames <- unique(sampleList[[1]])

holdingDir <- "1.platypusCalls/somaticTotal/"
namePrepended <- ".snv.annoVar.variant_function.txt"
newPrepended <- ".exc.snv.annoVar.variant_function.txt"
sharedPrep <- ".clonal.annoVar.variant_function.txt"

for(currSet in 1:length(setNames)){
  
  #subset main list
  subSample <- subset(sampleList, sampleList[1]==setNames[currSet])
  
  compName <- subSample[1,4]
  setName <- subSample[1,1]
  
  #get polyp vcf to filter
  polypvcf <- read.table(file=paste(subSample[1,6], holdingDir, setName,"/", setName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
  
  #get cancer vcf 
  canvcf <- read.table(file=paste(subSample[1,6], holdingDir, compName,"/", compName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
  
  polyRows <- c()
  polypCounter <- 1
  for(currRow in 1:nrow(polypvcf)){
    currChrom <- polypvcf[currRow, 3]
    currPos <- polypvcf[currRow, 4]
    
    #subset comparison vcf
    tempCompVcf <- canvcf[canvcf[[3]]==currChrom, ]
    
    if(currPos %in% tempCompVcf[[4]]){
      #variant is shared across cancer
    }else{
      #variant is polyp exclusive
      polyRows[polypCounter] <- currRow
      polypCounter <- polypCounter + 1
    }
  }
  
  #subset vcf and output
  filtPolypVcf <- polypvcf[polyRows, ]
  
  #get allele frequencies
  filtPolypVcf[12] <- filtPolypVcf[[10]] / filtPolypVcf[[8]]
  filtPolypVcf[13] <- filtPolypVcf[[11]] / filtPolypVcf[[9]]
  
  filtPolypVcf <- filtPolypVcf[filtPolypVcf[["V13"]] > 0.1, ]
  newVCFName <- paste(subSample[1,6], holdingDir, setName,"/", setName, newPrepended, sep="")
  write.table(filtPolypVcf, file=newVCFName, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
  #output shared variants to separate lists
  sharedPolypVcf <- polypvcf[-polyRows, ]
  
  #get allele frequencies
  sharedPolypVcf[12] <- sharedPolypVcf[[10]] / sharedPolypVcf[[8]]
  sharedPolypVcf[13] <- sharedPolypVcf[[11]] / sharedPolypVcf[[9]]
  
  sharedPolypVcf <- sharedPolypVcf[sharedPolypVcf[["V13"]] > 0.05, ]
  
  newVCFName <- paste(subSample[1,6], holdingDir, setName,"/", setName, sharedPrep, sep="")
  write.table(sharedPolypVcf, file=newVCFName, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  #count number of different elements and output table
  elemtTab <- table(sharedPolypVcf[[1]])
  newVCFName <- paste(subSample[1,6], holdingDir, setName,"/", setName, ".clonal.elements.txt", sep="")
  write.table(elemtTab, file=newVCFName, sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  
}


