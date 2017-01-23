# script calculates homoplasy vs non-homoplast differences
# plots ref depth vs non-ref depth
# extra: plots mutational burden boxplot


################### notes ###################



################### subroutines ####################



################### main program ####################

arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=1){
	stop("\usage: > Rscript diversityAnalysis.R < fileList.csv > < holding directory > < .prependedName.vcf > ")
}

#sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

#holdingDir <- arguments[2]
holdingDir <- "2.phylogenetics/exome.0.01/"
holdingDir2 <- "1.platypusCalls/exome.0.01/"

#vcfName <- arguments[3]
#varName <- ".snv.annoVar.variant_function.0.01.txt"
varName <- ".snv.annoVar.exonic_variant_function.0.01.txt"
homoplasyName <- ".homoAnno.txt"

#get samplenames list
sampleNames <- sampleList[c(1,2,7)]
sampleNamesList <- unique(sampleNames[1])

#remove unwanted columns and make unique set table
sampleList <- sampleList[-c(3:4)]
sampleList <- unique(sampleList[c(1,3:6)])

#list of binary tables for burden plot
listTabs <- as.list(NA)
totalSetNames <- unique(sampleList[[1]])

totalDivTab <- as.list(NA)
mainCounter <- 1
#change loop to 1:nrow(sampleList) for all samples
for(x in 1:length(totalSetNames)){
  #get variant file
  snvFileName <- paste(sampleList[x,3], holdingDir2, sampleList[x,1], "/",sampleList[x,1], varName, sep="")
  snvDataIn <- read.table(file=snvFileName, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
  snvDataIn[1] <- NULL
  
  #get homoplasy file
  homoFile <- paste(sampleList[x,3], holdingDir, sampleList[x,1], "/",sampleList[x,1], homoplasyName, sep="")
  homoData <- read.table(file=homoFile, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
  #homoData[1] <- paste(homoData[[3]], ":", homoData[[4]], sep="")
  
	#get sample names
	currentNames <- subset(sampleNames, sampleNames[1]==sampleList[x,1], select=V2)
	currentSet <- totalSetNames[x]
	
	noSamples <- sampleList[x,5]	
	normalIndex <- 1+sampleList[x,4]
  normalName <- currentNames[normalIndex, 1]
	
  #name columns
  names(snvDataIn) <- c("region", "gene", "chrom", "pos", "pos2", "alt", "ref", paste(currentNames[[1]], ".NR", sep=""), paste(currentNames[[1]], ".NV", sep=""))
  
	#remove normal columns
  snvDataIn[paste(normalName, ".NV", sep="")] <- NULL
  snvDataIn[paste(normalName, ".NR", sep="")] <- NULL
  currentNames <- currentNames[-normalIndex, 1]
	noSamples <- noSamples - 1
  
	#subtract total read depth from variant depth to get ref specific depth
	for(currCol in 1:noSamples){
	  snvDataIn[ncol(snvDataIn)+1] <- snvDataIn[[(currCol+7+noSamples)]] / snvDataIn[[(currCol+7)]]
    names(snvDataIn)[ncol(snvDataIn)] <- currentNames[currCol]
	}
  snvDataIn["pos2"] <- paste(snvDataIn[["chrom"]], ":", snvDataIn[["pos"]], sep="")
  
  #split main table into homoplasy and non-homoplasy
  homoConfData <- as.list(NA)
  homoConfData[[1]] <- snvDataIn[homoData[[1]],]
  homoConfData[[2]] <- snvDataIn[-homoData[[1]],] 
  names(homoConfData)[[1]] <- "homoplasy variants"
  names(homoConfData)[[2]] <- "non-homoplasy variants"
  
  #plot hostograms of VAF for homoplasy and non-homoplasy
  pdf(file=paste(sampleList[x,3], holdingDir, sampleList[x,1], "/",sampleList[x,1], "homoplasyVAFPlot.pdf", sep=""), height=(noSamples*4), width=10)
  par(mar=c(5,5,5,5), mfrow=c(noSamples, 2))
    for(currSam in 1:noSamples){
      for(currPlot in 1:2){
        plotData <- homoConfData[[currPlot]][[(currSam + 7 + (noSamples*2))]]
        plotData <- plotData[plotData > 0]
        hist(x=plotData, xlim=c(0,1), breaks = seq(0,1,0.01), xlab="VAF", main=paste(names(homoConfData)[[currPlot]], "sample", names(homoConfData[[currPlot]][(currSam + 7 + (noSamples*2))])) ) 
      }
    }
  dev.off()
  
  #plot histograms of depths for homoplasy and non-homoplasy
  pdf(file=paste(sampleList[x,3], holdingDir, sampleList[x,1], "/",sampleList[x,1], "homoplasyNRDepth.pdf", sep=""), height=(noSamples*4), width=10)
  par(mar=c(5,5,5,5), mfrow=c(noSamples, 2))
  for(currSam in 1:noSamples){
    for(currPlot in 1:2){
      hist(x=homoConfData[[currPlot]][[(currSam + 7)]], xlim=c(0,250), breaks = seq(0,250,5), xlab="NR read depth", main=paste(names(homoConfData)[[currPlot]], "sample", names(homoConfData[[currPlot]][(currSam + 7 + (noSamples*2))])) ) 
    }
  }
  dev.off()
  pdf(file=paste(sampleList[x,3], holdingDir, sampleList[x,1], "/",sampleList[x,1], "homoplasyNVDepth.pdf", sep=""), height=(noSamples*4), width=10)
  par(mar=c(5,5,5,5), mfrow=c(noSamples, 2))
  for(currSam in 1:noSamples){
    for(currPlot in 1:2){
      hist(x=homoConfData[[currPlot]][[(currSam + noSamples + 7)]], xlim=c(0,100), breaks = seq(0,2000,5), xlab="NV read depth", main=paste(names(homoConfData)[[currPlot]], "sample", names(homoConfData[[currPlot]][(currSam + 7 + (noSamples*2))])) ) 
    }
  }
  dev.off()
  
  
  #plot graphs of depth vs VAF
  pdf(file=paste(sampleList[x,3], holdingDir, sampleList[x,1], "/",sampleList[x,1], "homoplasyNRvsVAF.pdf", sep=""), height=(noSamples*4), width=10)
  par(mar=c(5,5,5,5), mfrow=c(noSamples, 2))
  for(currSam in 1:noSamples){
    VAFhomo <- homoConfData[[1]][[(currSam + (noSamples*2) + 7)]]
    VAFnonHomo <- homoConfData[[2]][[(currSam + (noSamples*2) + 7)]]
    NRhomo <- homoConfData[[1]][[(currSam + 7)]]
    NRnonHomo <- homoConfData[[2]][[(currSam + 7)]]  
    
    #remove zero values (variant not present in sample)
    NRhomo <- NRhomo[VAFhomo > 0]
    NRnonHomo <- NRnonHomo[VAFnonHomo > 0]
    VAFhomo <- VAFhomo[VAFhomo > 0]
    VAFnonHomo <- VAFnonHomo[VAFnonHomo > 0]
    
    #perform stats (null: samples are from the same population distribution)
    pvalVAF <- ks.test(x=VAFhomo, y=VAFnonHomo, alternative = "less")$p.value
    pvalNR <- ks.test(x=NRhomo, y=NRnonHomo, alternative = "less")$p.value
    
    for(currPlot in 1:2){
      yData <- homoConfData[[currPlot]][[(currSam + 7)]]
      xData <- homoConfData[[currPlot]][[(currSam + (noSamples*2) + 7)]]
      yData <- yData[xData > 0]
      xData <- xData[xData > 0]
      plot(x=xData, y=yData, pch=20, xlab="VAF", ylim=c(0,250), ylab="depth at loci", xlim=c(0,1), main=paste(names(homoConfData)[[currPlot]], "sample", names(homoConfData[[currPlot]][(currSam + 7 + (noSamples*2))])))
      if(currPlot == 2){
        text(x=0.8 , y=200, labels=paste("(ks.test, VAF) p =", round(pvalVAF, digit=3)))
        text(x=0.8 , y=180, labels=paste("(ks.test, depth) p =", round(pvalNR, digit=3)))
      }    
    }
  }
  dev.off()
}

