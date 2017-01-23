# script processes set.09 of the MSeq project exclusively 

########## info ###########


########## libraries ###########

#txtFileListSub <- rootNameFile; subListSub <- subList; holdingSub <- holdingDirVCF; outMutSigSub <- outSig
makeMutSig <- function(txtFileListSub, subListSub, holdingSub, outMutSigSub){
  #make mutSig table and output to .txt file
  mutSigTab <- data.frame(matrix(data=NA, nrow=1, ncol=4))
  names(mutSigTab) <- c("sample", "chrom", "pos", "mutation")
  
  for(y in txtFileListSub){
    #get vcf (.txt file) if contains variants
    fileTest <- file.info(paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", y ,sep=""))
    if(fileTest$size == 0){
      dataSet <- data.frame(matrix(0, nrow=1, ncol=4))
      names(dataSet) <- c("sample", "chrom", "pos", "mutation")
      dataSet[4] <- 0
      dataSet[3] <- 0
      dataSet[2] <- 0
      dataSet[1] <- y
      next
    }
    dataSet <- read.table(colClasses = "character", file=paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", y ,sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)
    
    #delete unwanted columns
    dataSet <- dataSet[,-c(5:ncol(dataSet))]
    dataSet[4] <- paste(dataSet[[3]], ">", dataSet[[4]], sep="")
    dataSet[3] <- dataSet[2]
    dataSet[2] <- as.character(dataSet[[1]])
    dataSet[1] <- y
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
    
    #add to mut signature table
    mutSigTab <- rbind(mutSigTab, dataSet)
  }
  
  #remove first row with NA
  mutSigTab <- mutSigTab[-1,]
  
  #remove non-conventional chromosomes and rename X and Y (if needed)
  mutSigTab[mutSigTab == "X"] <- 23
  mutSigTab[mutSigTab == "Y"] <- 24
  
  #remove non-number chromosomes
  mutSigTab <- subset(mutSigTab, chrom == 1 | chrom == 2 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10 | chrom == 11 | chrom == 12 | chrom == 13 | chrom == 14 | chrom == 15 | chrom == 16 | chrom == 17 | chrom == 18 | chrom == 19 | chrom == 20 | chrom == 21 | chrom == 22 | chrom == 23 | chrom == 24)
  
  #output mutSig table
  write.table(mutSigTab, file=paste(outMutSigSub, ".txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

#txtFileListSub <- txtFileList; subListSub <- subList; holdingSub<-holdingDirVCF; outMutSigSub <- outMutSig4
makeMutSigBanMerged <- function(txtFileListSub, subListSub, holdingSub, outMutSigSub){
  #make mutSig table and output to .txt file
  mutSigTab <- data.frame(matrix(data=NA, nrow=1, ncol=4))
  names(mutSigTab) <- c("sample", "chrom", "pos", "mutation")
  
  #subset txt file list for non-branch files
  branchTxtFileList <- txtFileListSub[grep("branch", txtFileList)]
  txtFileListSub <- txtFileListSub[-grep("branch", txtFileList)]
  
  for(y in txtFileListSub){
    fileTest <- file.info(paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", y, ".txt" ,sep=""))
    if(fileTest$size == 0){
      dataSet <- data.frame(matrix(0, nrow=1, ncol=4))
      names(dataSet) <- c("sample", "chrom", "pos", "mutation")
      dataSet[4] <- 0
      dataSet[3] <- 0
      dataSet[2] <- 0
      dataSet[1] <- y
      mutSigTab <- rbind(mutSigTab, dataSet)
      next
    }
    
    #get vcf (.txt file)
    dataSet <- read.table(colClasses = "character" ,file=paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", y, ".txt" ,sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)
    
    #delete unwanted columns
    dataSet <- dataSet[,-c(5:ncol(dataSet))]
    dataSet[4] <- paste(dataSet[[3]], ">", dataSet[[4]], sep="")
    dataSet[3] <- dataSet[2]
    dataSet[2] <- as.character(dataSet[[1]])
    dataSet[1] <- y
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
    
    #add to mut signature table
    mutSigTab <- rbind(mutSigTab, dataSet)
  }
  
  #remove first row with NA
  mutSigTab <- mutSigTab[-1,]
  
  #make temp mutSig table and output to .txt file
  mutSigBran <- data.frame(matrix(data=NA, nrow=1, ncol=4))
  names(mutSigBran) <- c("sample", "chrom", "pos", "mutation")
  
  #get merged branch files
  for(z in branchTxtFileList){
    fileTest <- file.info(paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", z, ".txt" ,sep=""))
    if(is.na(fileTest$size)){
      next
    }else if(fileTest$size == 0){
      next
    }
    
    #get vcf (.txt file)
    dataSet <- read.table(file=paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", z, ".txt" ,sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE)
    
    #delete unwanted columns
    dataSet <- dataSet[,-c(5:ncol(dataSet))]
    dataSet[4] <- paste(dataSet[[3]], ">", dataSet[[4]], sep="")
    dataSet[3] <- dataSet[2]
    dataSet[2] <- as.character(dataSet[[1]])
    dataSet[1] <- paste(subList[1,1], ".branches", sep="")
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
    
    #add to mut signature table
    mutSigBran <- rbind(mutSigBran, dataSet)
  }
  #remove duplicate entries in branch table
  mutSigBran <- mutSigBran[-1,]
  mutSigBran <- mutSigBran[order(mutSigBran[2]),]
  mutSigBran <- unique(mutSigBran)
  
  #merge tables
  if(nrow(mutSigBran) > 0){
    mutSigTab <- rbind(mutSigTab, mutSigBran)
  }
  
  #remove non-conventional chromosomes and rename X and Y (if needed)
  mutSigTab[mutSigTab == "X"] <- 23
  mutSigTab[mutSigTab == "Y"] <- 24
  
  #remove non-number chromosomes
  mutSigTab <- subset(mutSigTab, chrom == 1 | chrom == 2 | chrom == 2 | chrom == 3 | chrom == 4 | chrom == 5 | chrom == 6 | chrom == 7 | chrom == 8 | chrom == 9 | chrom == 10 | chrom == 11 | chrom == 12 | chrom == 13 | chrom == 14 | chrom == 15 | chrom == 16 | chrom == 17 | chrom == 18 | chrom == 19 | chrom == 20 | chrom == 21 | chrom == 22 | chrom == 23 | chrom == 24)
  
  #output mutSig table
  write.table(mutSigTab, file=paste(outMutSigSub, ".txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

plotSign <- function(subListSub, holdingSub, branName){
  #plotting variables
  cytoSig <- c("aCa", "aCc", "aCg", "aCt", "cCa", "cCc", "cCg", "cCt", "gCa", "gCc", "gCg", "gCt", "tCa", "tCc", "tCg", "tCt")
  thymineSig <- c("aTa", "aTc", "aTg", "aTt", "cTa", "cTc", "cTg", "cTt", "gTa", "gTc", "gTg", "gTt", "tTa", "tTc", "tTg", "tTt")
  mutSigVector <- c(cytoSig, cytoSig, cytoSig, thymineSig, thymineSig, thymineSig)
  
  mutMatrix <- paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", subListSub[1,1], branName, ".txt.mut.matrix", sep="")
  mutSamples <- paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", subListSub[1,1], branName, ".txt.samples", sep="")
  outMutAna <- paste(subListSub[1,6], holdingSub, subListSub[1,1], "/", subListSub[1,1], branName, ".txt", sep="")
  
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
  pdfName <- paste(subListSub[1,6], holdingDir, subListSub[1,1], "/", subListSub[1,1], branName, ".mutSig.raw.pdf", sep="")
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
  
  return(c(mutMatrix, mutSamples))
}

clusterData <- function(subListSub, holdingSub, branName, branSam){
  mutMatrix <- paste(subListSub[1,6], holdingSub, "/", branName, sep="")
  mutSamples <- paste(subListSub[1,6], holdingSub, "/", branSam, sep="")
  
  anaDataIn <- read.table(file=mutMatrix, header=FALSE, stringsAsFactors=FALSE)
  tumourNames <- read.table(file=mutSamples, header=FALSE, stringsAsFactors=FALSE)
  noSamples <- nrow(anaDataIn)
  row.names(anaDataIn) <- tumourNames[[1]]
  
  #for mutation matrix, produce proportions
  for(i in 1:nrow(anaDataIn)){
    tempSum <- sum(anaDataIn[i,])
    for(j in 1:ncol(anaDataIn)){
      anaDataIn[i,j] <- anaDataIn[i,j]/tempSum
    }
  }
  
  #invert the data.frame for bootstrap analysis
  anaDataInv <- t(anaDataIn)
  names(anaDataInv) <- row.names(anaDataIn)
  
  #perform clustering with bootstrap values (don't find same structure as below)
  dataClust <- pvclust(anaDataInv, method.hclust="ward.D", method.dist="euclidean", nboot=1000)
  pdfName <- paste(subListSub[1,6], holdingSub, "/", branName, ".mutSig.boot.pdf", sep="")
  pdf(file=pdfName, onefile=TRUE, width=12, height=8)
  plot(dataClust, cex=0.5)
  dev.off()
  
  #perfrom clustering, with default parameters
  distM <- dist(x=anaDataIn)
  dataClust <- hclust(d=distM)
  pdfName <- paste(subListSub[1,6], holdingSub, "/", branName, ".mutSig.clust.pdf", sep="")
  pdf(file=pdfName, onefile=TRUE, width=12, height=8)
  plot(dataClust, cex=0.5)
  dev.off()
  
  #perform clustering, with ward method and marked groupings
  dataClust <- hclust(d=distM, method="ward.D")
  pdfName <- paste(subListSub[1,6], holdingSub, "/", branName, ".mutSig.wardD.pdf", sep="")
  pdf(file=pdfName, onefile=TRUE, width=12, height=8)
  plot(dataClust, cex=0.5)
  groups <- cutree(dataClust, k=8)
  rect.hclust(dataClust, k=8, border="red")
  dev.off()	
}

mergeMutSig <- function(mutFileListSub, outNameSub){
  #merge sample name and matrix files
  mutFiles <- paste(mutFileListSub, collapse=" ")
  totalOutFile <- paste(subList[1,6], holdingDir, outNameSub, sep="")	
  system(command=paste("cat ", mutFiles, " > ", totalOutFile, sep=""))
}



#### main program ####


#requires connection to apocrita

subList <- read.csv(file="~/PhD/CRCproject/sampleList.Set.09.filtered.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- "1.platypusCalls/somaticTotal/"
namePrepended <- ".snv.annoVar.variant_function"

holdingDirVCF <- "1.platypusCalls/somaticTotal/"

noInfoCol <- 7

#perform pipeline for thi set only
normalIndex <- (subList[1,7]+1)
normalName <- subList[normalIndex,2]
samplesLocs <- c(1:subList[1,8])
samplesLocs <- samplesLocs[-normalIndex]
noSamples <- subList[1,8]

#get somatic conformance file
annoName <- paste(subList[1,6], holdingDir, subList[1,1], "/", subList[1,1], namePrepended, ".txt", sep="")
confData <- read.table(file= annoName, sep="\t", header=FALSE, stringsAsFactors=FALSE)
confData <- confData[confData[[3]]=="1" | confData[[3]]=="2" | confData[[3]]=="3" | confData[[3]]=="4" | confData[[3]]=="5" | confData[[3]]=="6" | confData[[3]]=="7" | confData[[3]]=="8" | confData[[3]]=="9" | confData[[3]]=="10" | confData[[3]]=="11" | confData[[3]]=="12" | confData[[3]]=="13" | confData[[3]]=="14" | confData[[3]]=="15" | confData[[3]]=="16" | confData[[3]]=="17" | confData[[3]]=="18" | confData[[3]]=="19" | confData[[3]]=="20" | confData[[3]]=="21" | confData[[3]]=="22" | confData[[3]]=="X" | confData[[3]]=="Y", ]
names(confData) <- c("type", "genes", "chrom", "pos", "pos2", "ref", "alt", paste(subList[[2]], ".NR", sep=""), subList[[2]])

confData["Set9_84"] <- 0


#define branch and leafs to be measured
branchList <- data.frame(matrix(NA, nrow = 21, ncol=14))
names(branchList) <- c("ID", "phylogenetic region",	"number of variants",	"downstream", "exonic",	"intergenic",	"intronic",	"ncRNA_exonic",	"ncRNA_intronic",	"ncRNA_UTR3",	"splicing",	"upstream",	"UTR3",	"UTR5")
rownames(branchList) <- c("root", "distal carcinoma and both intervening", "distal carcinoma and distal intervening", "distal intervening", "proximal intervening", "distal carcinoma root", "distal carcinoma B1", "distal carcinoma B2", "distal biopsy 1", "distal biopsy 2", "distal biopsy 3", "distal biopsy 4", "distal biopsy 5", "proximal carcinoma root",  "proximal carcinoma B1", "proximal carcinoma B2", "proximal biopsy 1", "proximal biopsy 2", "proximal biopsy 3", "proximal biopsy 4", "proximal biopsy 5")
branchList[1, 2] <- paste(c(1:6, 8:13), collapse = ":")
branchList[2, 2] <- paste(c(1:6, 13), collapse = ":")
branchList[3, 2] <- paste(c(1:6), collapse = ":")
branchList[4, 2] <- paste(c(6), collapse = ":")
branchList[5, 2] <- paste(c(13), collapse = ":")
branchList[6, 2] <- paste(c(1:5), collapse = ":")
branchList[7, 2] <- paste(c(2:5), collapse = ":")
branchList[8, 2] <- paste(c(3:5), collapse = ":")
branchList[9, 2] <- paste(c(1), collapse = ":")
branchList[10, 2] <- paste(c(2), collapse = ":")
branchList[11, 2] <- paste(c(3), collapse = ":")
branchList[12, 2] <- paste(c(4), collapse = ":")
branchList[13, 2] <- paste(c(5), collapse = ":")
branchList[14, 2] <- paste(c(8:12), collapse = ":")
branchList[15, 2] <- paste(c(9:12), collapse = ":")
branchList[16, 2] <- paste(c(10,12), collapse = ":")
branchList[17, 2] <- paste(c(8), collapse = ":")
branchList[18, 2] <- paste(c(9), collapse = ":")
branchList[19, 2] <- paste(c(10), collapse = ":")
branchList[20, 2] <- paste(c(11), collapse = ":")
branchList[21, 2] <- paste(c(12), collapse = ":")

branchList["ID"] <- c("root", "D_inter", "D_Dinter", "Dinter", "Pinter", "DcarRoot", "DcarB1", "DcarB2", "Dcar1", "Dcar2", "Dcar3", "Dcar4", "Dcar5", "PcarRoot",  "PcarB1", "PcarB2", "Pcar1", "Pcar2", "Pcar3", "Pcar4", "Pcar5")

dataColumns <- c(1:13) + noInfoCol + noSamples

#populate table and output mutational signatures
for(currRow in 1:nrow(branchList)){
  print(paste("##### assessing branch:", branchList[currRow, "ID"], "#######"))
  
  testVect <- rep(0,13)
  currCols <- as.numeric(strsplit(branchList[currRow, "phylogenetic region"], split = ":")[[1]])
  testVect[currCols] <- 1
  
  assessRows <- c() 
  for(currKeep in 1:nrow(confData)){
    assTab <- confData[currKeep, dataColumns] > 0
    assTab[assTab == "TRUE"] <- 1
    assTab[assTab == "FALSE"] <- 0
    assTab <- as.vector(assTab)
    
    if(identical(assTab, testVect)){
      assessRows[currKeep] <- "TRUE"
    }else{
      assessRows[currKeep] <- "FALSE"
    }
  }
  filtConf <- confData[assessRows=="TRUE",]
  
  #populate table
  branchList[currRow, "number of variants"] <- nrow(filtConf)
  tabFiltConf <- table(filtConf[["type"]])
  branchList[currRow, "downstream"] <- as.numeric(tabFiltConf["downstream"])
  branchList[currRow, "exonic"] <- as.numeric(tabFiltConf["exonic"])
  branchList[currRow, "intergenic"] <- as.numeric(tabFiltConf["intergenic"])
  branchList[currRow, "intronic"] <- as.numeric(tabFiltConf["intronic"])
  branchList[currRow, "ncRNA_exonic"] <- as.numeric(tabFiltConf["ncRNA_exonic"])
  branchList[currRow, "ncRNA_intronic"] <- as.numeric(tabFiltConf["ncRNA_intronic"])
  branchList[currRow, "ncRNA_UTR3"] <- as.numeric(tabFiltConf["ncRNA_UTR3"])
  branchList[currRow, "splicing"] <- as.numeric(tabFiltConf["splicing"])
  branchList[currRow, "upstream"] <- as.numeric(tabFiltConf["upstream"])
  branchList[currRow, "UTR3"] <- as.numeric(tabFiltConf["UTR3"])
  branchList[currRow, "UTR5"] <- as.numeric(tabFiltConf["UTR5"])

  filteredVCFTxt <- paste(subList[1,6], holdingDir, subList[1,1], "/", subList[1,1], ".", branchList[currRow, "ID"], ".txt", sep="")
  write.table(filtConf, file=filteredVCFTxt, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  #remove unwanted columns
  filtConf <- filtConf[-c(1:2,5)]
  filtTxtInput <- paste(subList[1,6], holdingDir, subList[1,1], "/", subList[1,1], ".", branchList[currRow, "ID"], ".Emu.txt", sep="")
  write.table(filtConf, file=filtTxtInput, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  #make mutational signature
  rootNameFile <- paste(subList[1,1], ".", branchList[currRow, "ID"], ".Emu.txt", sep="")
  outSig <- paste(subList[1,6], holdingDir, subList[1,1], "/", subList[1,1], ".", branchList[currRow, "ID"], ".sig", sep="")
  makeMutSig(rootNameFile, subList, holdingDirVCF, outSig)
  system(command=paste("~/bin/EMu-prepare --chr /Users/cross01/mount/william/referenceHG19/chromosomeFastas/ --mut ", outSig,".txt" ,sep=""), wait=TRUE)
  
  #plot graph of signature
  mutMatFiles <- plotSign(subList, holdingDir, paste(".", branchList[currRow, "ID"], ".sig", sep=""))
  
}


#output varaints table
branchTableOut <- "~/Dropbox/MSeq-CRC/mseq.writeup/suppTables/supp.table.07-synchronousVariants.csv"
write.csv(branchList, file = branchTableOut, quote = FALSE, row.names = TRUE, col.names = TRUE)


