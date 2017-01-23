# plot summary of CNV disruptions from MSI cancer

##########################   notes   ##########################

########################## libraries  ##########################

######################### subroutines ##########################

#filter valid CNV events and assign phylogenetic location ####
#segDataSub <- segCNV;noSamSub <- noSamples
filterCNVs <- function(segDataSub, noSamSub){
  returnList <- as.list(NA)
  
  #assess segmentations if normal, remove
  removeSegs <- c()
  counterSub <- 1
  segDataSub <- as.data.frame(append(segDataSub, x = NA))
  names(segDataSub)[1] <- "phyloLoc"
  for(currSeg in 1:nrow(segDataSub)){
    assessSeg <- table(segDataSub[currSeg, 6:(5+(noSamSub*2))] == 1)
    if("TRUE" %in% names(assessSeg)){
      if(as.integer(assessSeg["TRUE"])==(2*noSamSub)){
        #don't keep, its diploid
        removeSegs[counterSub] <- currSeg
        counterSub <- counterSub + 1
      }else{
        #keep
      }
    }
  }
  #subset to remove normal regions
  if(length(removeSegs)>0){
    segDataSub <- segDataSub[-removeSegs,]
  }
  
  #assess phylogenetic location
  for(currSeg in 1:nrow(segDataSub)){
    #get current CNV seg into table
    assessVector <- segDataSub[currSeg, 6:(5+(noSamSub*2))]
    assessTable <- data.frame(matrix(NA, ncol=2, nrow=noSamSub, byrow = TRUE))
    seqAp <- seq(1, (noSamSub*2), 2)
    for(currAp in 1:nrow(assessTable)){
      assessTable[currAp, 1] <- assessVector[seqAp[currAp]]
      assessTable[currAp, 2] <- assessVector[seqAp[currAp]+1]
    }
    assessTable[3] <- paste(assessTable[[1]], ":", assessTable[[2]], sep="")
    names(assessTable) <- c("minor", "major", "merged")
    
    #assess phylo by uniqueness of states
    stateTest <- table(assessTable[["merged"]])
    noStates <- length(unique(names(stateTest)))
    
    #assess number of states
    if(noStates == 1){
      #trunkal
      segDataSub[currSeg, 1] <- "T"
    }else if("1:1" %in% names(stateTest)){
      if(stateTest["1:1"] == (noSamSub-1)){
        #true leaf event
        segDataSub[currSeg, 1] <- "L"
      }else if(noStates == 2){
        #its a true branch event
        segDataSub[currSeg, 1] <- "B"
      }else{
        #its a bifuricated branch (two branch events in same region)
        segDataSub[currSeg, 1] <- "BiB"
      }
    }else{
      #if no samples have WT state they are ether branch or trunk
      if(noStates == 2){
        #its a bifuricated trunk (one branch events in same region as a trunk)
        segDataSub[currSeg, 1] <- "BiT"
      }else{
        #phylogenetic loc cannot be determined
        segDataSub[currSeg, 1] <- "unknown >2 bifurcations"
      }
    }
  }
  returnList[[1]] <- segDataSub
  returnList[[2]] <- removeSegs
  return(returnList)
}

#function returns cytoband annotation
#subTab <- annoTable[currAnno,]; subBandList <- cytobandList; subDrivers <- annoDriverList
getCytoAnno <- function(subTab, subBandList, subDrivers){
  #subset drivers list
  subsettedDrive <- subDrivers[subDrivers[[1]]==as.numeric(subTab[1, "chr"]), ]
  startEnd <- subTab[1, c("first.locus", "last.locus")]
  
  #subset cyto bands file
  subsetBands <- subBandList[subBandList[[1]]==as.numeric(subTab[1, "chr"]), ]
  subsetBands <- subsetBands[subsetBands[[3]]>as.numeric(startEnd["first.locus"]) & subsetBands[[3]]<=as.numeric(startEnd["last.locus"]), ]
  
  if(nrow(subsetBands)==1){
    #get possible driver genes
    driverGenesSub <- subsettedDrive[subsettedDrive[[2]]>=subsetBands[1,2] & subsettedDrive[[3]]<=subsetBands[1,3], ]
    
    annotation <- subsetBands[1, 4] 
  }else{
    #get possible driver genes
    driverGenesSub <- subsettedDrive[subsettedDrive[[2]]>=subsetBands[1,2] & subsettedDrive[[3]]<=subsetBands[nrow(subsetBands),3], ]
    
    annotation <- paste(subsetBands[1, 4], ":", subsetBands[nrow(subsetBands), 4], sep="") 
  }
  annotationFinal <- as.list(NA)
  annotationFinal[[1]] <- annotation
  if(nrow(driverGenesSub)>0){
    annotationFinal[[2]] <- driverGenesSub
  }else{
    annotationFinal[[2]] <- NA
  }
  
  return(annotationFinal)
}

#### get SNV for each region ####
getSNVs <- function(currPosStartSub, currPosEndSub, currChromSub, currSamSub, vcfDataSub){
  #subset by chromosome location
  vcfDataSub <- subset(vcfDataSub, vcfDataSub$chrom==currChromSub)
  keepListSub <- c()
  counterSub <- 1
  
  #get current sample AF col name
  currSampCol <- paste(currSamSub, ".AF", sep="")
  
  #loop though variants and keep those in CNV
  for(currRow in 1:nrow(vcfDataSub)){
    if(vcfDataSub[currRow, "pos"] >= currPosStartSub & vcfDataSub[currRow, "pos"] <= currPosEndSub & vcfDataSub[currRow, currSampCol] > 0){
      #keep variant as it is in CNV region
      keepListSub[counterSub] <- currRow
      counterSub <- counterSub + 1
    }
  }
  vcfDataSub <- vcfDataSub[keepListSub, c("chrom", "pos", "ref", "alt", currSampCol)]
  return(vcfDataSub)
}

#dataSet <- clonalSNVs; IDname <- annoTable[currRow, "ID"]; holdingSub <- holdingDir; outMutSigSub <- outSig
makeMutSig <- function(dataSet, IDname, holdingSub, outMutSigSub){
 
  #delete unwanted columns
  dataSet <- dataSet[,-c(9:ncol(dataSet))]
  dataSet[1] <- IDname
  dataSet[2] <- dataSet[4]
  dataSet[3] <- dataSet[5]
  dataSet[4] <- paste(dataSet[[7]], ">", dataSet[[8]], sep="")
  dataSet <- dataSet[1:4]
  names(dataSet) <- c("sample", "chrom", "pos", "mutation")
  
  #correct for chr names if needed
  chromList <- strsplit(dataSet[[2]], "chr")
  if(length(chromList[[1]]) == 1){
    #do nothing
  }else{
    #change names from chrX to X 
    for(w in 1:length(chromList)){
      dataSet[w,2] <- chromList[[w]][2]
    }		
  }
  
  #remove non-conventional chromosomes and rename X and Y (if needed)
  dataSet[dataSet == "X"] <- 23
  dataSet[dataSet == "Y"] <- 24
  
  #remove non-number chromosomes
  dataSet <- subset(dataSet, chrom == 1 | chrom == 2 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10 | chrom == 11 | chrom == 12 | chrom == 13 | chrom == 14 | chrom == 15 | chrom == 16 | chrom == 17 | chrom == 18 | chrom == 19 | chrom == 20 | chrom == 21 | chrom == 22 | chrom == 23 | chrom == 24)
  
  #output mutSig table
  write.table(dataSet, file=paste(outMutSigSub, ".txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

#subListSub <- subSample; holdingSub <- holdingDir; branName <- outSig
plotSign <- function(subListSub, holdingSub, branName){
  #plotting variables
  cytoSig <- c("aCa", "aCc", "aCg", "aCt", "cCa", "cCc", "cCg", "cCt", "gCa", "gCc", "gCg", "gCt", "tCa", "tCc", "tCg", "tCt")
  thymineSig <- c("aTa", "aTc", "aTg", "aTt", "cTa", "cTc", "cTg", "cTt", "gTa", "gTc", "gTg", "gTt", "tTa", "tTc", "tTg", "tTt")
  mutSigVector <- c(cytoSig, cytoSig, cytoSig, thymineSig, thymineSig, thymineSig)
  
  mutMatrix <- paste(branName, ".txt.mut.matrix", sep="")
  mutSamples <- paste(branName, ".txt.samples", sep="")
  outMutAna <- paste(branName, ".txt", sep="")
  
  #plot raw mutational signatures
  anaDataIn <- read.table(file=mutMatrix, header=FALSE, stringsAsFactors=FALSE)
  tumourNames <- read.table(file=mutSamples, header=FALSE, stringsAsFactors=FALSE)
  noSamples <- nrow(anaDataIn)
  
  #convert anaDataIn to matrix
  temp2 <- as.matrix(anaDataIn)
  temp <- matrix(NA, nrow=96, ncol=noSamples)
  for(j in 1:noSamples){
    temp[,j] <- temp2[j,]
  }
  
  noSamples <- nrow(tumourNames)
  if(noSamples <=12 ){
    plotHeight <- 30
  }
  if(noSamples >12){
    plotHeight <- 50
  }
  if(noSamples >40){
    plotHeight <- 80
  }
  
  #pdf file name
  pdfName <- paste(branName, ".mutSig.raw.pdf", sep="")
  pdf(file=pdfName, onefile=TRUE, width=12, height=plotHeight)
  if(noSamples <13){
    par(mfrow=c(12,1), xpd=TRUE, mar=c(3.5,2.5,2.5,2.5))
  }
  if(noSamples>12){
    par(mfrow=c(24,1), xpd=TRUE, mar=c(3.5,2.5,2.5,2.5))
  }
  if(noSamples>24){
    par(mfrow=c(24,2), xpd=TRUE, mar=c(3.5,2.5,2.5,2.5))
  }
  if(noSamples>40){
    par(mfrow=c(48,2), xpd=TRUE, mar=c(3.5,2.5,2.5,2.5))
  }
  
  for(j in 1:noSamples){
    barplot(temp[,j], xaxt="n", col=c(rep("blue",16), rep("black",16), rep("red",16), rep("grey",16), rep("green",16), rep("pink",16)), main=tumourNames[j,1])
    axis(side=1, at=seq(0.5,114.5,length.out=96), labels=mutSigVector, las=2, lty=0, cex.axis=0.5)
  }
  
  #add custom legend
  plot.new()
  legend('bottom',legend=c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), col=c("blue","black","red", "grey", "green", "pink"), lty=1, xpd=TRUE, cex = 0.75, lwd = 3)
  dev.off()
}


######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

setNames <- "Set.04"

#working dir
CNVholdingDir <- "7.CNVcalls.final/baf.updated/"
CNVfileNames <- ".penalty0.95.baf.gt.txt"

#SNV dir
variantCallFile <- ".snv.annoVar.variant_function.0.01.txt"
SNVworkingDir <- "1.platypusCalls/somaticTotal.0.01/"

#out directory
CNVout <- "7.CNVcalls.final/CNVannotation/"

#annotation lists
cytobandList <- read.table(file="~/PhD/ReferenceGenome/cytoBand.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")
annoDriverList <- read.csv(file="~/PhD/CRCproject/9.driverEvents/archive/driverList.chromLoci.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- "5.mutationalSignatures/"

############ 1.make segmentations summary table ###############

for(currSam in 1:length(setNames)){
  #subset main list
  subSample <- subset(sampleList, sampleList[1]==setNames[currSam])
  
  sampleNames <- subSample[[2]]
  sampleNames <- sampleNames[-(subSample[1,7]+1)]
  normalSample <- subSample[1,7]+1
  noSamples <- length(sampleNames)
  
  #get seg file
  segCNV <- read.table(file=paste(subSample[1,6], CNVholdingDir, subSample[1,1], CNVfileNames, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
  
  #set ambigous calls to WT
  segCNV[segCNV=="X"] <- 1
  
  #save sampled genome length then filter small regions
  segCNV["nloci"] <- (segCNV[["last.locus"]] - segCNV[["first.locus"]]) / 1000000

  #ensure states are numeric
  for(currChg in 5:ncol(segCNV)){
    segCNV[currChg] <- as.numeric(segCNV[[currChg]])
  }
  
  returnedData <- filterCNVs(segCNV, noSamples)
  diploidySegs <- segCNV[returnedData[[2]], ]
  diploidChromosomes <- unique(diploidySegs[[1]])
  filteredSegCNVs <- returnedData[[1]]

  #copy data to get annotation table
  annoTable <- filteredSegCNVs
  
  #annotate and rearragen table
  sortOrder <- c("T", "BiT", "B", "BiB", "L", "unknown >2 bifurcations")
  #annoTable <- annoTable[annoTable[["nloci"]] >= 5,]
  annoTable <- annoTable[order(match(annoTable[[1]], sortOrder)), ]
  annoTable["cytoBand"] <- NA
  annoTable["driverGenes"] <- NA
  
  #add cytoband
  for(currAnno in 1:nrow(annoTable)){
    annoReturn <- getCytoAnno(annoTable[currAnno,], cytobandList, annoDriverList)
    annoTable[currAnno, "cytoBand"] <- annoReturn[[1]][1]
    
    if(length(annoReturn[[2]]) != 0){
      annoTable[currAnno, "driverGenes"] <- paste(unlist(annoReturn[[2]][5]), collapse = ":")
    }else{
      annoTable[currAnno, "driverGenes"] <- NA
    }
  }
  
  #output annotation table
  #outAnno <- paste(subSample[1,6], CNVout, subSample[1,1], ".CNV.anno.txt", sep="")
  #write.table(annoTable, file=outAnno, sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    
  vcfFileName <- paste(subSample[1,6], SNVworkingDir, subSample[1,1], "/", subSample[1,1], variantCallFile, sep="")
  vcfData <- read.table(file=vcfFileName, sep="\t", stringsAsFactors = FALSE, header=FALSE)
  
  #remove normal sample
  vcfData <- vcfData[-c((8 + normalSample) , (8 + (noSamples + 1) + normalSample))]
  
  #calculate allele frequencies and change X/Y to 23/24
  for(y in 1:noSamples){
    vcfData[(ncol(vcfData)+1)] <- as.vector(vcfData[(8+(noSamples)+y)]) / as.vector(vcfData[8+y])
  }
  names(vcfData) <- c("line", "region", "variant", "chrom", "pos", "pos2", "ref", "alt", paste(sampleNames, ".NR", sep=""), paste(sampleNames, ".NV", sep=""), paste(sampleNames, ".AF", sep=""))
  chromData <- vcfData[["chrom"]]
  chromData[chromData=="X"] <- 23
  chromData[chromData=="Y"] <- 24
  vcfData["chrom"] <- chromData
  
  #get blank table 
  excTab <- data.frame(0, nrow=0, ncol=3)

  vcfDataUnfilt <- vcfData
  
  #remove snv if not truncal or in exclusion table above
  remRow <- c()
  remCounter <- 1
  for(currSNV in 1:nrow(vcfData)){
    SNVchrom <- vcfData[currSNV, 4]
    SNVpos <- vcfData[currSNV, 5]
    
    #if no CNV region exists for this chromosome, remove
    subCNV <- filteredSegCNVs[filteredSegCNVs[["chr"]]==SNVchrom, ]
    if(nrow(subCNV) == 0){
      remRow[remCounter] <- currSNV
      remCounter <- remCounter + 1
      next()
    }
    
    #if there are CNV regions to filter, test current variant
    subHyper <- excTab[excTab[[1]]==SNVchrom, ] 
    if(nrow(subHyper) != 0){
      #if truncal, see if variant is in filter chromosome regions (subHyper)
      for(currHyp in 1:nrow(subHyper)){
        if(SNVpos >= subHyper[currHyp, 2] & SNVpos <= subHyper[currHyp, 3]){
          remRow[remCounter] <- currSNV
          remCounter <- remCounter + 1
          next()
        }
      }
    }
    
    #variant not truncal, remove
    if(as.numeric(table(vcfData[currSNV, paste(sampleNames, ".AF", sep="")] > 0)["TRUE"]) != noSamples){
      remRow[remCounter] <- currSNV
      remCounter <- remCounter + 1
      next()
    }
    
    #finally, see if SNV in a CNV region
    for(currTest in 1:nrow(subCNV)){
      if(SNVpos >= subCNV[currTest, "first.locus"] & SNVpos <= subCNV[currTest, "last.locus"]){
        break()
      }else if(currTest == nrow(subCNV)){
        #remove as not in any CNV
        remRow[remCounter] <- currSNV
        remCounter <- remCounter + 1
      }
    }
  }
  
  vcfNew <- vcfData[-remRow, ]
  
  
  #get allele frequency for a diploid region of the genome
  remRowDip <- c()
  remDipCounter <- 1
  for(currSNV in 1:nrow(vcfData)){
    SNVchrom <- vcfData[currSNV, 4]
    SNVpos <- vcfData[currSNV, 5]
    
    #if no CNV region exists for this chromosome, remove
    subCNV <- diploidySegs[diploidySegs[["chr"]]==SNVchrom, ]
    if(nrow(subCNV) == 0){
      remRow[remCounter] <- currSNV
      remCounter <- remCounter + 1
      next()
    }
    
    #variant not truncal, remove
    if(as.numeric(table(vcfData[currSNV, paste(sampleNames, ".AF", sep="")] > 0)["TRUE"]) != noSamples){
      remRowDip[remDipCounter] <- currSNV
      remDipCounter <- remDipCounter + 1
      next()
    }
    
    #finally, see if SNV in a CNV region
    for(currTest in 1:nrow(subCNV)){
      if(SNVpos >= subCNV[currTest, "first.locus"] & SNVpos <= subCNV[currTest, "last.locus"]){
        break()
      }else if(currTest == nrow(subCNV)){
        #remove as not in any CNV
        remRowDip[remDipCounter] <- currSNV
        remDipCounter <- remDipCounter + 1
      }
    }
  }
  vcfDiploid <- vcfData[-remRowDip, ]
  
  #add annotation columns to table
  annoTable["noClonalVar"] <- NA
  annoTable["ID"] <- NA
  annoTable["ID"] <- paste(annoTable[["chr"]], "-", annoTable[["first.locus"]], sep="")
  
  removeCols <- paste(sampleNames, "_Major", sep="")
  annoTable <- annoTable[,-which(names(annoTable) %in% removeCols)]
  
  #this loop assumed there is only one SCNA per chromosome
  for(currRow in 1:nrow(annoTable)){
    #clonal or subclonal?
    if(annoTable[currRow, "phyloLoc"] == "T"){
      #clonal varaints in this SCNA region
      clonalSNVs <- vcfNew[vcfNew[["chrom"]]==annoTable[currRow, "chr"], ]
      annoTable[currRow, "noClonalVar"] <- nrow(clonalSNVs)
      
      #get mutational signature
      outSig <- paste(subSample[1,6], holdingDir, subSample[1,1], "/", subSample[1,1], ".", annoTable[currRow, "ID"], ".sig", sep="")
      makeMutSig(clonalSNVs, annoTable[currRow, "ID"], holdingDir, outSig)
      
      #use EMu to get character bins
      system(command=paste("~/bin/EMu-prepare --chr /Users/cross01/mount/william/referenceHG19/chromosomeFastas/ --mut ", outSig,".txt" ,sep=""), wait=TRUE)
      
      #plot graph of signature
      plotSign(subSample, holdingDir, outSig)
      
      #plot VAFs for this regions against diploid regions
      VAFplot <- paste(outSig, ".VAF.pdf", sep="")
      pdf(file=VAFplot, width = 10, height = (noSamples*5))
      par(mfrow=c(noSamples,2))
        for(currBio in 1:length(sampleNames)){
          currBioName <- sampleNames[currBio]
          VAFColTemp <- paste(currBioName, ".AF", sep="") 
          currVAFs <- clonalSNVs[[VAFColTemp]]
          comparisonVAFs <- vcfDiploid[[VAFColTemp]]
          hist(currVAFs, breaks = seq(0,1,0.05), xlab = "variant allele frequency", main=currBioName, cex.axis = 1.5)
          abline(v = mean(currVAFs), col="red")
          hist(comparisonVAFs, breaks = seq(0,1,0.05), xlab = "variant allele frequency", main=paste(currBioName, "diploid regions"), cex.axis = 1.5)
          abline(v = mean(comparisonVAFs), col="red")
          
          #perform stats
          statsRes <- wilcox.test(currVAFs, comparisonVAFs)
          legend("topright", legend = paste("wilcon; ", statsRes$p.value))
          
          annoTable[currRow, paste(currBioName, "_Minor", sep="")] <- statsRes$p.value
        }
      dev.off()
      
    }else{
      #subclonal, do nothing
    }
  }
  
  #output annotation table
  outTableName <- paste(subSample[1,6], holdingDir, subSample[1,1], "/", subSample[1,1], ".CNVs.txt", sep="")
  write.table(annoTable, file=outTableName, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}
