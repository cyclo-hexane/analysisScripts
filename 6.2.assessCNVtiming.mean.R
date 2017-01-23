# takes modelled CNV timing and plots figures
# adds estimates of truncal timing
# this version plots means of the sample timings

################### notes ###################
#
#   6.0.CNVtvaraintClustering.R (est alpha/beta and make input table: X.matlab.input.txt)
#             |
#             V
#   AnalyseCNVData.m (matlab script on apocrita, output: X.matlab.output.txt)
#             |
#             V
#   6.1.estimateBifurcations.R (used for plot)
#             |
#             V
#   6.2.assessCNVtiming.heatmap.R (make final figure of timings)
#
################# subroutines #################

#function returns cytoband annotation
#subTab <- timingInfo[1,]; subBandList <- cytobandList; subDrivers <- annoDriverList; startEnd <- subParam
getCytoAnno <- function(subTab, subBandList, subDrivers, startEnd){
  #subset drivers list
  subsettedDrive <- subDrivers[subDrivers[[1]]==as.numeric(subTab[1, "Chromosome"]), ]
  
  #subset cyto bands file
  subsetBands <- subBandList[subBandList[[1]]==as.numeric(subTab[1, "Chromosome"]), ]
  subsetBands <- subsetBands[subsetBands[[2]]>=as.numeric(startEnd["start"]) & subsetBands[[3]]<=as.numeric(startEnd["end"]), ]

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

############### main program ################

#libraries
library(plotrix)

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

#timing sample names
#sampleNames <- c("Set.02", "Set.03", "Set.05", "Set.09.Distal", "Set.09.Proximal", "Set.10", "Polyp.08.WGS")
#sampleNames <- "Set.10.recalled"
sampleNames <- c("Set.03", "Set.05", "Set.09.Distal", "Set.09.Proximal", "Set.10", "Polyp.08.WGS")

#load driver gene list
annoDriverList <- read.csv(file="~/PhD/CRCproject/9.driverEvents/archive/driverList.chromLoci.csv", header=FALSE, stringsAsFactors=FALSE)

#biomartGeneList <- read.table(file="~/PhD/ReferenceGenome/mart_export_genes.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
cytobandList <- read.table(file="~/PhD/ReferenceGenome/cytoBand.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")

CNVtimingDir <- "6.CNVtiming/modelResults.final/"
timingPrep <- ".matlab.output.txt"

#input data (only used to get start-end regions for CNVs)
SNVdataDir <- "6.CNVtiming/inputData.final/"
SNVApp <- ".input.parameters.ammend.txt"

#get bifurcation estimates
#biFile <- "~/PhD/CRCproject/6.CNVtiming/bifurcationPoints/bifurcationEstimates.csv"
#bifurcationTable <- read.csv(file=biFile, stringsAsFactors = FALSE, header=TRUE, row.names = 1)

#for each sample set  make timing graph by chromosome
for(currCan in 1:length(sampleNames)){
  
  currentSample <- subset(sampleList, sampleList[1]==sampleNames[currCan])
  nonNormSam <- currentSample[[2]]
  nonNormSam <- nonNormSam[-(currentSample[1,7]+1)]
  noSamples <- length(nonNormSam)
  
  #get CNV and timings
  timingFile <- paste(currentSample[1,6], CNVtimingDir, currentSample[1,1], "/", currentSample[1,1], timingPrep, sep="")
  CNVtotal <- read.table(file=timingFile, header=TRUE, stringsAsFactors = FALSE, sep=" ")
  CNVtotal["length"] <- CNVtotal[["length"]] / 1000000
  
  #get parameters file
  paramFile <- paste(currentSample[1,6], SNVdataDir, currentSample[1,1], "/", currentSample[1,1], SNVApp, sep="")
  paramData <- read.table(file=paramFile, header=TRUE, stringsAsFactors = FALSE, sep=",")
  paramData["CNVlength"] <- paramData[["CNVlength"]] / 1000000
  
  #get plotting parameters
  noPlotRows <- 23
  timingInterval <- c(0, 0.000004)

  
  ######### plot graph of total timing #########
  
  #colour list, tetrasomy, trisomy, LOH
  colList <- c("red4", "red", "blue")
  
  outputLoc <- paste(currentSample[1,6], CNVtimingDir, currentSample[1,1], "/", currentSample[1,1], ".CNVtiming.mean.pdf", sep="")
  pdf(file=outputLoc, width=12, height=5)
  par(xpd = TRUE, cex.main = 2.5, xpd=TRUE, mar=c(6,2,2,2)) 
  
 
  ########### main timing plot ###########
  plot(1, 1, col="white", axes=F, ylim=c(0, 24), xlim=timingInterval, xlab="", ylab="", main="")
  
  #counters
  plotCounter <- 0
  
  for(currChr in 1:23){    

    #add chromosome name
    text(x = rep(-0.00000012, 2), y = rep(plotCounter, 2), labels = paste("chr", currChr), cex=1)
    
    #add sample trace line
    lines(x = c(0, 0.000004), y=rep(plotCounter, 2), lwd=0.25, lty=2, col="grey35")
    
    #get samples with timings for this chromosome
    currSub <- CNVtotal[CNVtotal[["Chromosome"]]==currChr, ]
    
    if(nrow(currSub) == 0){
      plotCounter <- plotCounter +1
      next
    }
    
    #no CNVs to plot
    noCNVs <- unique(currSub[["length"]])
    
    #annotation table
    annoTab <- data.frame(matrix(NA, nrow=length(noCNVs), ncol=4))
    colnames(annoTab) <- c("CNV_ID", "anno", "pos", "col")

    #annotation table counter
    annoCounter <- 1
    
    tempPlotCounter <- plotCounter
    
    #plot each timing value and confidence interval
    for(currPoint in length(noCNVs):1){
      
      #subset parameter list
      subParam <- unique(paramData[paramData[["CNVlength"]]==noCNVs[currPoint], c("start", "end")])
      
      #subset chromosome list for this CNV ID and get mean values
      timingInfo <- currSub[currSub[["length"]]==noCNVs[currPoint], ]
      timingInfo["t1_est_single"] <- mean(timingInfo[["t1_est_joint"]])
      timingInfo["t1_MSE_single"] <- mean(timingInfo[["t1_MSE_joint"]])
      
      #add copy state column
      timingInfo["copyState"] <- paste(timingInfo[["a2"]], ":", timingInfo[["b2"]], sep="")

      #get cytoband and effected driver genes
      annotationInfo <- getCytoAnno(timingInfo[1,], cytobandList, annoDriverList, subParam)
      splitAnno <- strsplit(annotationInfo[[1]][1], split = ":")
      if(length(splitAnno[[1]]) > 1){
        staLoci<- substr(splitAnno[[1]][1], start = 1, stop = 1)
        endLoci<- substr(splitAnno[[1]][2], start = 1, stop = 1)
      }else{
        staLoci<- substr(splitAnno[[1]][1], start = 1, stop = 1)
        endLoci<- substr(splitAnno[[1]][1], start = 1, stop = 1)
      }
      
      #determine plot shape by q or p arm status
      if(staLoci == "p" & endLoci == "q"){
        #if both plot as triangle
        bioRep <- 17
      }else if(staLoci == "p" & endLoci == "p"){
        #if p arm plot as square
        bioRep <- 15
      }else{
        #else mark circle as q arm
        bioRep <- 18
      }
      
      
      #save information to annotation table
      annoTab[annoCounter, 1] <- timingInfo[1, "CNV_ID"]
      annoTab[annoCounter, 2] <- annotationInfo[[1]]
      annoTab[annoCounter, 3] <- timingInfo[1, "t1_est_joint"]
        
      #assign colour by CNV type
      if(timingInfo[1, "copyState"] == "2:2"){
        tempCol <- colList[1]
      }else if(timingInfo[1, "copyState"] == "1:2"){
        tempCol <- colList[2]
      }else{
        tempCol <- colList[3]
      }
      
      #add colour to anno tab
      annoTab[annoCounter, 4] <- tempCol
      annoCounter <- annoCounter + 1
      
      #mean confidence interval
      upperInter <- timingInfo[1, "t1_est_single"] + (sqrt(timingInfo[1, "t1_MSE_single"]) * 2)
      lines(x=c(timingInfo[1, "t1_est_single"], upperInter), y=rep(tempPlotCounter, 2), col=tempCol, lwd=0.8, lty=1)
      
      lowerInter <- timingInfo[1, "t1_est_single"] - (sqrt(timingInfo[1, "t1_MSE_single"]) * 2)
      lines(x=c(timingInfo[1, "t1_est_single"], lowerInter), y=rep(tempPlotCounter, 2), col=tempCol, lwd=0.8, lty=1)
      
      #main timing point
      points(x=timingInfo[1, "t1_est_single"], y=tempPlotCounter, pch=bioRep, lty = 0, col=tempCol, cex=2.3)
      
      tempPlotCounter <- tempPlotCounter + 0.1
    }
    
    
    #intrement counter at add annotation space
    plotCounter <- plotCounter + 1
    
    #add cytoband annotations
    CNVlist <- unique(annoTab[[1]])
    annoCNVlist <- data.frame(matrix(NA, nrow=0, ncol=3))
    for(currCNV in 1:length(CNVlist)){
     tempTab <- annoTab[annoTab[[1]]==CNVlist[currCNV], ]
     annoCNVlist <- rbind(tempTab[nrow(tempTab),], annoCNVlist)
    }
    text(x=annoCNVlist[[3]], y=rep((plotCounter-3), nrow(annoCNVlist)), labels = annoCNVlist[[2]], cex=0.3, srt=45, pos=3, col="black")
     
  }
  plotFrom <- plotCounter
  plotCounter <- plotCounter + 5
  
  #add sampling time
  sampleingIntervals <- unique(CNVtotal[["T_est_joint"]])
  samplingInt <- mean(sampleingIntervals)
  lines(x=rep(samplingInt, 2), y=c(0, 23), col="grey45", lwd=3, lty=1)
  
  #add x-axis
  axis(1, at=seq(0, 0.000004, 0.0000005), labels=seq(0, 0.000004, 0.0000005), lwd=1, las=2, cex.axis=1, line=0.5)
  
  dev.off()
  
  #convert timing table to mean values
  tempCNVID <- unlist(strsplit(CNVtotal[["CNV_ID"]], split = "-"))
  numSeq <- c(1:length(tempCNVID))
  numSeq <- numSeq[-seq(1,length(numSeq),4)]
  tempCNVID <- tempCNVID[numSeq]
  counter <- 1
  for(i in seq(1, length(tempCNVID), 3)){
    CNVtotal[counter, "CNV_ID"] <- paste(tempCNVID[i], "-", tempCNVID[i+1], "-", tempCNVID[i+2], sep="") 
    counter <- counter + 1
  }
  Idlist <- unique(CNVtotal[["CNV_ID"]])
  
  t1Means <- c()
  for(currID in 1:length(Idlist)){
    t1Means[currID] <- mean(CNVtotal[CNVtotal[["CNV_ID"]]==Idlist[currID], "t1_est_joint"])
  }
  
  #make QQplots for timings
  outputLoc <- paste(currentSample[1,6], CNVtimingDir, currentSample[1,1], "/", currentSample[1,1],".qqplots.CNVtiming.pdf", sep="")
  pdf(file=outputLoc, width=10, height=5)
  par(mfrow=c(1, 2), mar=c(3,3,2,2))
  
  qqnorm(t1Means, main = paste("Normal QQ", currentSample[1,1]))
  qqline(t1Means, col="red")
  
  hist(t1Means, main ="histogram of timings", breaks=seq(0, 0.000004, 0.000000025))
  
  #test for normality, null= came from a normal dist
  #par(xpd=TRUE)
  pValue <- shapiro.test(t1Means)$p.value
  text(y=2, x=0.0000035, labels = paste("p.value = ", round(pValue, digits = 5), sep=""), cex=0.75)

  dev.off()
  
  
  
}

