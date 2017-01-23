# makes evolution schematics
# takes input sample list and filtering (removal) list as arguments
# uses, exome annovar, cloneHD segmentations, and driver list data
# 1.determines CNV driver distruptions and trunk proportions
# 2.determines SNV driver distributions
# 3.plots evolution schematic
 
################### notes ###################

#biomartGeneList has been subsetted to only contain driver genes (making searches quicker)

################# libraries #################


############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#remSampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.total.exclusions.csv", header=FALSE, stringsAsFactors=FALSE)
SNVholdingDir <- "1.platypusCalls/totalData/"
SNVPrepended <- ".total.filt.0.05.chr.txt"
CNVholdingDir <- "7.CNVcalls.final/baf.updated/"
CNVPrepended <- ".phyloCNVs.csv"

sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[-c(21:24)]

#load driver gene list
#annoDriverList <- read.csv(file="~/PhD/CRCproject/9.driverGenes/driverList.chromLoci.csv", header=FALSE, stringsAsFactors=FALSE)

#biomartGeneList <- read.table(file="~/PhD/ReferenceGenome/mart_export_genes.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
cytobandList <- read.table(file="~/PhD/ReferenceGenome/cytoBand.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")

#annotate each set each sample 
for(currSamp in 1:length(sampleNames)){
	print(paste("#### making schematic for sample ", sampleNames[currSamp], " ####",sep=""))
	
	#subset main list and get info
	subSample <- subset(sampleList, sampleList[1]==sampleNames[currSamp])
	#remSubSample <- subset(remSampleList, remSampleList[1]==sampleNames[currSamp])
  setName <- unique(subSample[[1]])
	noSamples <- subSample[1,8]-1
	normalIndex <- 1+subSample[1,7]

  
	################## 1.CNV driver analysis ################
  
	#setup input/output names and get files
	CNVfile <- paste(subSample[1,6], CNVholdingDir, setName, CNVPrepended, sep="")
	CNVdataIn <- read.csv(file=CNVfile, header=TRUE, stringsAsFactors=FALSE, colClasses = c("character"))
  
  #get driver genes within CNV regions
	CNVdataIn[ncol(CNVdataIn) + 1] <- NA
  names(CNVdataIn)[ncol(CNVdataIn)] <- "genes"
	CNVdataIn[ncol(CNVdataIn) + 1] <- NA
	names(CNVdataIn)[ncol(CNVdataIn)] <- "loci"
  
  #loop through and populate CNVdataIn with gene names
  for(currCNV in 1:nrow(CNVdataIn)){
    currChrom <- CNVdataIn[currCNV, "chr"]
    startCNV <- as.numeric(CNVdataIn[currCNV, "first.locus"])
    endCNV <- as.numeric(CNVdataIn[currCNV, "last.locus"])
    
    #get banding loci for current interval
    lociList <- c()
    lociCounter <- 1
    cytoBandSub <- cytobandList[cytobandList[1]==currChrom,]
    if(nrow(cytoBandSub)==0){
      CNVdataIn[currCNV, ncol(CNVdataIn)-1] <- "none"
      #CNVdataIn[currCNV, ncol(CNVdataIn)] <- "none"
      next
    }
    cytoAnno <- paste(cytoBandSub[1,1], cytoBandSub[1,4], ":", cytoBandSub[nrow(cytoBandSub),4], sep="")
    annoDriverSub <- annoDriverList[annoDriverList[1]==currChrom, ]
    if(nrow(annoDriverSub)!=0){
      for(currLoci in 1:nrow(cytoBandSub)){
        if(cytoBandSub[currLoci,2] > startCNV & cytoBandSub[currLoci,3] < endCNV){
          lociList[lociCounter] <- as.character(cytoBandSub[currLoci,4])
          lociCounter <- lociCounter + 1
        }
      }
      
      #get gene list using position
      tempGeneList <- c()
      geneCounter <- 1
      for(currLoci in 1:nrow(annoDriverSub)){
        if(annoDriverSub[currLoci, 2] >= startCNV & annoDriverSub[currLoci, 3] <= endCNV){
          tempGeneList[geneCounter] <- currLoci
          geneCounter <- geneCounter + 1
        }
      }
      CNVdataIn[currCNV, ncol(CNVdataIn)] <- cytoAnno
    }else{
      tempGeneList <- c()
      CNVdataIn[currCNV, ncol(CNVdataIn)] <- "none"
    }
    
    #add results to main table (CNVdataIn)
    if(!is.null(tempGeneList)){
      CNVdataIn[currCNV, ncol(CNVdataIn)-1] <- paste(annoDriverSub[tempGeneList, 5], collapse = ":")
    }else{
      CNVdataIn[currCNV, ncol(CNVdataIn)-1] <- "none"
    }
  }
	
  
  #sort into trunk - branch - leaf
	trunkRows <- which(CNVdataIn["phyloLoc"] == "T")
	branchRows <- which(CNVdataIn["phyloLoc"] == "B")
	leafRows <- which(CNVdataIn["phyloLoc"] == "L")
	orderRows <- c(trunkRows, branchRows, leafRows)
	CNVdataIn <- CNVdataIn[orderRows,]
  
  #output annotated files
  fileNameOut <- paste(subSample[1,6], CNVholdingDir, setName, ".phyloCNVs.anno.csv", sep="")
  write.table(CNVdataIn, file=fileNameOut, sep=",", quote = FALSE, row.names=FALSE, col.names = FALSE)
}  
  
