# takes modelled CNV timing and plots figures
# adds estimates of truncal timing
# this version plots means of the sample timings and stats of linearity across time

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
sampleNames <- c("Set.03", "Set.05", "Set.09.Distal", "Set.09.Proximal", "Set.10", "Polyp.08.WGS")

#load driver gene list
annoDriverList <- read.csv(file="~/PhD/CRCproject/9.driverEvents/archive/driverList.chromLoci.csv", header=FALSE, stringsAsFactors=FALSE)

#biomartGeneList <- read.table(file="~/PhD/ReferenceGenome/mart_export_genes.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
cytobandList <- read.table(file="~/PhD/ReferenceGenome/cytoBand.txt", header=FALSE, stringsAsFactors=FALSE, sep="\t")

CNVtimingDir <- "6.CNVtiming/modelResults.final/"
timingPrep <- ".matlab.output.txt"

timingTabList <- as.list(NA)

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

  CNVtotal["CNV_ID"] <- paste(CNVtotal[["Chromosome"]], "-", CNVtotal[["length"]], sep="")
  uniqueIDs <- unique(CNVtotal[["CNV_ID"]])
  
  #get mean timings for each ID
  newTimingTab <- data.frame(matrix(NA, nrow=length(uniqueIDs), ncol=3))
  newTimingTab[1] <- uniqueIDs
  for(x in 1:length(uniqueIDs)){
    newTimingTab[x, 2] <- mean(CNVtotal[CNVtotal[["CNV_ID"]]==uniqueIDs[x], "t1_est_joint"])
    newTimingTab[x, 3] <- unique(CNVtotal[CNVtotal[["CNV_ID"]]==uniqueIDs[x], "length"])
  }
  timingTabList[[currCan]] <- newTimingTab[order(newTimingTab[[2]]), ]
  names(timingTabList)[[currCan]] <- currentSample[1,1]
  
  #setup breaks table
  breaksInt <- seq(0, 0.000004, 0.00000005)
  plotBreaksTab <- data.frame(matrix(NA, nrow=length(breaksInt), ncol=3))
  plotBreaksTab[1] <- breaksInt
  for(currCl in 1:nrow(plotBreaksTab)){
    plotBreaksTab[currCl, 2] <- nrow(newTimingTab[newTimingTab[[2]] < plotBreaksTab[currCl, 1], ])
    if(plotBreaksTab[currCl, 2] != 0){
      plotBreaksTab[currCl, 3] <- plotBreaksTab[currCl, 2] / length(uniqueIDs) 
    }else{
      plotBreaksTab[currCl, 3] <- 0
    }
  }
  
  #get stats 
  tempVect <- newTimingTab[[2]]
  print(ks.test(seq(from = min(tempVect), to = max(tempVect), length.out = length(tempVect)), tempVect))
  plot(seq(from = min(tempVect), to = max(tempVect))
  points()
  
} 
  
  ######### plot cumulative graph of total timing #########
  
  
  outputLoc <- paste(currentSample[1,6], CNVtimingDir, currentSample[1,1], "/", currentSample[1,1], ".CNVtiming.cumulative.pdf", sep="")
  pdf(file=outputLoc, width=12, height=5)
  par(xpd = TRUE, cex.main = 2.5, xpd=TRUE, mar=c(6,2,2,2)) 
    plot(x=plotBreaksTab[[1]], y=plotBreaksTab[[2]], pch=20, cex=0.5)
    lines(plotBreaksTab[[1]], plotBreaksTab[[2]], col="red")
    rect(xleft=CNVtotal[1, "T_est_joint"], xright=CNVtotal[1, "T_est_joint"]+0.000000005, ybottom = 0, ytop = 15)
    text(0.0000001, 1, labels = paste("p=", signif(timingStats$p.value, digits = 2)))
  dev.off()
  
}


#plot against time with linear regression
outputLoc <- paste("~/Dropbox/MSeq-CRC/mseq.writeup/suppFigures/supp.figure.06-CNVsizeVtiming.pdf", sep="")
pdf(file=outputLoc, width=10, height=15)
par(xpd = FALSE, cex.main = 2.5, mar=c(4,4,4,4), mfrow=c(3,2)) 
for(currSam in 1:length(sampleNames)){
  currTab <- timingTabList[[sampleNames[currSam]]]
  lmResults <- lm(data = currTab, formula = X3 ~ X2)
  nd <- data.frame(x=seq(0, (2*max(currTab[["X2"]])), length=100))
  names(nd) <- "X2"
  pConf1 <- predict(lmResults, interval="prediction", newdata=nd)
  plot(x=currTab[[2]], y=currTab[[3]], pch=20, cex=1.5, cex.axis=1.5, xlab = "evolutionary time", ylab = "size (Mb)", main=sampleNames[currSam])
  abline(lmResults, col="red")
  matlines(nd, pConf1[,c("lwr","upr")], col="grey50", lty=2, type="p", pch=1, cex=0.2)
  text(0.000002, 30, labels = paste("p=", round(summary(lmResults)$coefficients[2,4], digits = 4)))
  text(0.000002, 28, labels = paste("R2=", round(summary(lmResults)$adj.r.squared, digits = 4)))
}
dev.off()



#plot size against time with linear regression for all merged
outputLoc <- paste("~/Dropbox/MSeq-CRC/mseq.writeup/suppFigures/supp.figure.06-CNVsizeVtiming.b.pdf", sep="")
pdf(file=outputLoc, width=10, height=15)
par(xpd = FALSE, cex.main = 2.5, mar=c(4,4,4,4), mfrow=c(3,2)) 


  currTab <- rbind(timingTabList[[1]], timingTabList[[2]], timingTabList[[3]], timingTabList[[4]], timingTabList[[5]], timingTabList[[6]])
  lmResults <- lm(data = currTab, formula = X3 ~ X2)
  nd <- data.frame(x=seq(0, (2*max(currTab[["X2"]])), length=100))
  names(nd) <- "X2"
  pConf1 <- predict(lmResults, interval="prediction", newdata=nd)
  plot(x=currTab[[2]], y=currTab[[3]], pch=20, cex=1.5, cex.axis=1.5, xlab = "evolutionary time", ylab = "size (Mb)", main=sampleNames[currSam])
  abline(lmResults, col="red")
  matlines(nd, pConf1[,c("lwr","upr")], col="grey50", lty=2, type="p", pch=1, cex=0.2)
  text(0.000002, 30, labels = paste("p=", round(summary(lmResults)$coefficients[2,4], digits = 4)))
  text(0.000002, 28, labels = paste("R2=", round(summary(lmResults)$adj.r.squared, digits = 4)))

dev.off()


