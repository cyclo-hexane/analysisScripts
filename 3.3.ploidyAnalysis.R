# statistically assesses trunk and non-trunk CNV size
# additionally assesses distrubution of loss, gain CNVs

##########################   notes   ##########################

########################## libraries  ##########################

######################### subroutines ##########################

######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
sampleList <- sampleList[-c(125:142), ]

sampleNames <- unique(sampleList[[1]])
noSets <- length(sampleNames)

normalCols <- c()
normalCounter <- 1
for(currSet in 1:length(sampleNames)){
  currTab <- sampleList[sampleList[[1]]==sampleNames[currSet], ]
  tempNorCol <- which(currTab[(1 + currTab[1,7]), 2] == sampleList[[2]])
  for(currNo in 1:length(tempNorCol)){
    normalCols[normalCounter] <- tempNorCol[currNo] 
    normalCounter <- normalCounter + 1 
  }
}
normalCols <- unique(normalCols)
sampleList <- sampleList[-normalCols, ]

#working dir
CNVholdingDir <- "7.CNVcalls.archive/baf.updated/"
CNVfileNames <- ".penalty0.95.baf.gt.txt"
plotDir <- "10.finalFigures/"


##################### acquire CNV mutational counts ##########################
CNVTable <- data.frame(matrix(NA, ncol=9, nrow=0))
names(CNVTable) <- c("Set", "biopsy", "diploid", "triploid", "tetraploid", "loss", "Bloss", "polyploid", "av.ploidy")

for(j in 97:nrow(sampleList)){
	print(paste("#### geting CNV info for sample ", sampleList[j, 2], " ####",sep=""))
  
  tempTab <- data.frame(matrix(NA, ncol=9, nrow=1))
  names(tempTab) <- c("Set", "biopsy", "diploid", "triploid", "tetraploid", "loss", "Bloss", "polyploid", "av.ploidy")
  
  setName <- sampleList[j,1]
  
  #get normal biopsy name for this set
  subList <- sampleList[sampleList[[1]]==setName, ]
	currBiopsy <- sampleList[j,2]
	
	tempTab[1, 1] <- setName
	tempTab[1, 2] <- currBiopsy
	
	if(j > 96){
	  currMaj <- paste("X",currBiopsy, "_Major", sep="")
	  currMin <- paste("X",currBiopsy, "_Minor", sep="")
	}else{
	  currMaj <- paste(currBiopsy, "_Major", sep="")
	  currMin <- paste(currBiopsy, "_Minor", sep="")
	}
	
	#setup input/output names
	dataIn <- read.table(file=paste(sampleList[j,6], CNVholdingDir, setName, CNVfileNames, sep=""), header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep="\t")
	
	#subset to current biopsy
	dataIn <- dataIn[, c("chr", "first.locus", "nloci", "last.locus", currMin, currMaj)]
	dataIn[dataIn == "X"] <- 1
	
	dataIn[3] <- dataIn[[4]] - dataIn[[2]]
	totalRegions <- sum(as.numeric(dataIn[[3]]))
	
	#get average ploidy for this sample based on a 10kb window
	total10kbRegions <- 0
	averagePloidy <- 0
	for(currCNV in 2:nrow(dataIn)){
	  tempRegionInver <- dataIn[currCNV, "nloci"] / 10000
	  total10kbRegions <- total10kbRegions + tempRegionInver
	  averagePloidy <- averagePloidy + (tempRegionInver * (as.numeric(dataIn[currCNV, 5]) + as.numeric(dataIn[currCNV, 6])))
	}
	resPloidy <- averagePloidy / total10kbRegions
	tempTab[1, "av.ploidy"] <- resPloidy
	
	
	#loop through segmentations and determine ploidy
	lossTotal <- 0
	BlossTotal <- 0
	tetraTotal <- 0
	triTotal <- 0
	dipTotal <- 0
	polyTotal <- 0
	
	for(currRow in 1:nrow(dataIn)){
	  if(dataIn[currRow, 5] == 1 & dataIn[currRow, 6] == 1){
	    #diploid region
	    dipTotal <- dipTotal + dataIn[currRow, 3]
	  }else if(dataIn[currRow, 5] == 0 & dataIn[currRow, 6] == 1){
	    #loss region
	    lossTotal <- lossTotal + dataIn[currRow, 3]
	  }else if(dataIn[currRow, 5] == 0 & dataIn[currRow, 6] > 1){
	    #Bloss region
	    BlossTotal <- BlossTotal + dataIn[currRow, 3]
	  }else if(dataIn[currRow, 6] > 3){
	    polyTotal <- polyTotal + dataIn[currRow, 3]
	  }else if(dataIn[currRow, 5] == 1 & dataIn[currRow, 6] == 2){
	    #triploid region
	    triTotal <- triTotal + dataIn[currRow, 3]
	  }else{
	    #tetraploid region
	    tetraTotal <- tetraTotal + dataIn[currRow, 3]
	  }
	}
	
	#add to final table
	tempTab[1, "diploid"] <- dipTotal / totalRegions
	tempTab[1, "triploid"] <- triTotal / totalRegions
	tempTab[1, "tetraploid"] <- tetraTotal / totalRegions
	tempTab[1, "loss"] <- lossTotal / totalRegions
	tempTab[1, "Bloss"] <- BlossTotal / totalRegions
	tempTab[1, "polyploid"] <- polyTotal / totalRegions
	CNVTable <- rbind(CNVTable, tempTab)

	#if this is last sample in set get add mean ploidy and values to CNVTableMean
	CNVtemp <- CNVTable[CNVTable[[1]]==setName, ]
	tempTabMean <- data.frame(matrix(NA, ncol=9, nrow=1))
	names(tempTabMean) <- c("Set", "biopsy", "diploid", "triploid", "tetraploid", "loss", "Bloss", "polyploid", "av.ploidy")
	tempTabMean[1,1] <- setName
	tempTabMean[1,2] <- "mean"
	if(j == nrow(sampleList)){
	  tempTabMean[1, "diploid"] <- mean(CNVtemp[["diploid"]])
	  tempTabMean[1, "triploid"] <- mean(CNVtemp[["triploid"]])
	  tempTabMean[1, "tetraploid"] <- mean(CNVtemp[["tetraploid"]])
	  tempTabMean[1, "loss"] <- mean(CNVtemp[["loss"]])
	  tempTabMean[1, "Bloss"] <- mean(CNVtemp[["Bloss"]])
	  tempTabMean[1, "polyploid"] <- mean(CNVtemp[["polyploid"]])
	  tempTabMean[1, "av.ploidy"] <- mean(CNVtemp[["av.ploidy"]])
	  CNVTable <- rbind(CNVTable, tempTabMean)
	}else if(setName != sampleList[(j+1), 1]){
	  tempTabMean[1, "diploid"] <- mean(CNVtemp[["diploid"]])
	  tempTabMean[1, "triploid"] <- mean(CNVtemp[["triploid"]])
	  tempTabMean[1, "tetraploid"] <- mean(CNVtemp[["tetraploid"]])
	  tempTabMean[1, "loss"] <- mean(CNVtemp[["loss"]])
	  tempTabMean[1, "Bloss"] <- mean(CNVtemp[["Bloss"]])
	  tempTabMean[1, "polyploid"] <- mean(CNVtemp[["polyploid"]])
	  tempTabMean[1, "av.ploidy"] <- mean(CNVtemp[["av.ploidy"]])
	  CNVTable <- rbind(CNVTable, tempTabMean)
	}
	
}


########### plot means only with error bars ########### 


#get mean values and order table
CNVMeans <- CNVTable[CNVTable[[2]]=="mean", ]
sortVector <- c(rep("C", 11), rep("A", 9))
CNVMeans <- CNVMeans[order(sortVector, CNVMeans[[3]]), ]
CNVMeans <- CNVMeans[c("Set", "biopsy", "diploid", "polyploid", "tetraploid", "triploid", "Bloss", "loss", "av.ploidy")]

pdf(file=paste(sampleList[1,6], plotDir, "ploidyPlot.mean.pdf", sep=""), onefile=TRUE, width=8, height=5)
par(xpd=TRUE, mar=c(4,6,2,2))

#start plotting
plot(1, 1, col="white", axes=F, xlim=c(0,6.6), ylim=c(0, length(sampleNames)), xlab="", ylab="", main="")

#colour list
colList <- c("grey", "purple", "red3", "red", "blue", "lightsteelblue")

#for each sample set plot rectangles
plotCounter <- 1
for(currPlot in 1:nrow(CNVMeans)){
  colTemp <- colList
  #dip rect
  rect(xleft = 0, xright = CNVMeans[plotCounter, 3], ybottom = plotCounter, ytop = plotCounter+1, col=colTemp[1], border = "white")
  
  #poly rect
  rect(xleft = 1.1, xright = 1.1+CNVMeans[plotCounter, 4], ybottom = plotCounter, ytop = plotCounter+1, col=colTemp[2], border = "white")
  
  #tetra rect
  rect(xleft = 2.2, xright = 2.2+CNVMeans[plotCounter, 5], ybottom = plotCounter, ytop = plotCounter+1, col=colTemp[3], border = "white")
  
  #tri rect
  rect(xleft = 3.3, xright = 3.3+CNVMeans[plotCounter, 6], ybottom = plotCounter, ytop = plotCounter+1, col=colTemp[4], border = "white")
  
  #LOH rect
  rect(xleft = 4.4, xright = 4.4+CNVMeans[plotCounter, 7], ybottom = plotCounter, ytop = plotCounter+1, col=colTemp[5], border = "white")
  
  #LOH rect
  rect(xleft = 5.5, xright = 5.5+CNVMeans[plotCounter, 8], ybottom = plotCounter, ytop = plotCounter+1, col=colTemp[6], border = "white")
  
  
  #dd average ploidy
  text(x=6.6, y=(plotCounter+0.5), label=round(CNVMeans[plotCounter, "av.ploidy"], 3), pos = 4, cex = 0.5)
  
  
  text(x=round(CNVMeans[plotCounter, 3], 3), y=(plotCounter+0.5), label=round(CNVMeans[plotCounter, 3], 3), pos = 4, cex = 0.3)
  text(x=(1.1 + round(CNVMeans[plotCounter, 4], 3)), y=(plotCounter+0.5), label=round(CNVMeans[plotCounter, 4], 3), pos = 4, cex = 0.3)
  text(x=(2.2 + round(CNVMeans[plotCounter, 5], 3)), y=(plotCounter+0.5), label=round(CNVMeans[plotCounter, 5], 3), pos = 4, cex = 0.3)
  text(x=(3.3 + round(CNVMeans[plotCounter, 6], 3)), y=(plotCounter+0.5), label=round(CNVMeans[plotCounter, 6], 3), pos = 4, cex = 0.3)
  text(x=(4.4 + round(CNVMeans[plotCounter, 7], 3)), y=(plotCounter+0.5), label=round(CNVMeans[plotCounter, 7], 3), pos = 4, cex = 0.3)
  text(x=(5.5 + round(CNVMeans[plotCounter, 8], 3)), y=(plotCounter+0.5), label=round(CNVMeans[plotCounter, 8], 3), pos = 4, cex = 0.3)
  
  
  #add confidence lines
  currName <- CNVMeans[plotCounter, 1]
  currConf <- CNVTable[CNVTable[[1]]==currName, ]
  tempDip <- currConf[["diploid"]]
  tempTet <- currConf[["tetraploid"]]
  tempTri <- currConf[["triploid"]]
  tempLOH <- currConf[["Bloss"]]
  tempLoss <- currConf[["loss"]]
  temppoly <- currConf[["polyploid"]]
  
  #add confidence lines to graph
  lines(x=c(min(tempDip), max(tempDip)), y=c(plotCounter+0.5, plotCounter+0.5), col="black", lwd=0.8)
  lines(x=c(1.1+min(temppoly), 1.1+max(temppoly)), y=c(plotCounter+0.5, plotCounter+0.5), col="black", lwd=0.8)
  lines(x=c(2.2+min(tempTet), 2.2+max(tempTet)), y=c(plotCounter+0.5, plotCounter+0.5), col="black", lwd=0.8)
  lines(x=c(3.3+min(tempTri), 3.3+max(tempTri)), y=c(plotCounter+0.5, plotCounter+0.5), col="black", lwd=0.8)
  lines(x=c(4.4+min(tempLOH), 4.4+max(tempLOH)), y=c(plotCounter+0.5, plotCounter+0.5), col="black", lwd=0.8)
  lines(x=c(5.5+min(tempLoss), 5.5+max(tempLoss)), y=c(plotCounter+0.5, plotCounter+0.5), col="black", lwd=0.8)
  
  plotCounter <- plotCounter + 1
}  

#add proportion axis
axis(side=1, at=seq(0,1,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.4, las=2, line=-1.5)
axis(side=1, at=seq(1.1,2.1,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.4, las=2, line=-1.5)
axis(side=1, at=seq(2.2,3.2,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.4, las=2, line=-1.5)
axis(side=1, at=seq(3.3,4.3,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.4, las=2, line=-1.5)
axis(side=1, at=seq(4.4,5.4,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.4, las=2, line=-1.5)

#add names
axis(side=2, at=seq(1.5,(length(sampleNames)+0.5), 1), labels=CNVMeans[[1]], cex.axis=0.4, las=2, tick = FALSE, line=-1.5)

#add ploidy names
axis(side=3, at=seq(0.5, 6.6, 1), labels=c("diploid", "polyploid", "tetraploid", "triploid", "LOH", "loss", "av.ploidy"), cex.axis=0.8, tick = FALSE, line=0)

#add lines
lines(x=c(1.1,1.1),  y=c(1, length(sampleNames)+1), col="black")
lines(x=c(2.2,2.2),  y=c(1, length(sampleNames)+1), col="black")
lines(x=c(3.3,3.3),  y=c(1, length(sampleNames)+1), col="black")
lines(x=c(4.4,4.4),  y=c(1, length(sampleNames)+1), col="black")
lines(x=c(5.5,5.5),  y=c(1, length(sampleNames)+1), col="black")
lines(x=c(0,0),  y=c(1, length(sampleNames)+1), col="black")

dev.off()



############# plot tetraploidy correlations ###############

CNVTet <- CNVTable[CNVTable[[2]]!="mean", ]
tetCNVs <- CNVTet[CNVTet[[1]]=="Set.10" | CNVTet[[1]]=="Set.09.Distal" | CNVTet[[1]]=="Set.03" | CNVTet[[1]]=="Set.08", ] 
nonTetraCNVs <- CNVTable[CNVTable[[2]]!="mean" & CNVTable[[1]]!="Set.10" & CNVTable[[1]]!="Set.09.Distal" & CNVTable[[1]]!="Set.03" & CNVTable[[1]]!="Set.08", ]
nonTetraCNVs <- nonTetraCNVs[c(1:52), ]

pdf(file=paste(sampleList[1,6], plotDir, "tetraploidCorr.pdf", sep=""), onefile=TRUE, width=8, height=4)
par(xpd=FALSE, mar=c(6,6,2,2), mfrow=c(1,2))
  #plot data the produce correlations
  
  #tetraploidy vs trisomy etc in suspected tetraploidy sets
  plot(main="carcinomas with < 10% diploid", tetCNVs$tetraploid, tetCNVs$triploid, xlim = c(0,1), ylim=c(0,1), pch=20, col="red", xlab = "proportion of genome at tetraploidy", ylab = "proportion of genome")
  points(tetCNVs$tetraploid, tetCNVs$Bloss, col="blue", pch=20)
  points(tetCNVs$tetraploid, tetCNVs$loss, col="lightsteelblue", pch=20)
  
  #fir regression lines
  R2fit1 <- lm(triploid ~ tetraploid, data=tetCNVs)
  abline(R2fit1, col="red")
  R2fit2 <- lm(Bloss ~ tetraploid, data=tetCNVs)
  abline(R2fit2, col="blue")
  R2fit3 <- lm(loss ~ tetraploid, data=tetCNVs)
  abline(R2fit3, col="lightsteelblue")
  
  text(x=0.6, y=0.9, paste("R2:", format(summary(R2fit1)$r.squared, digits=4), "p=", round(summary(R2fit1)$coefficients[2,4], 12)), col="red")
  text(x=0.6, y=0.8, paste("R2:", format(summary(R2fit2)$r.squared, digits=4), "p=", round(summary(R2fit2)$coefficients[2,4], 12)), col="blue")
  text(x=0.6, y=0.7, paste("R2:", format(summary(R2fit3)$r.squared, digits=4), "p=", round(summary(R2fit3)$coefficients[2,4], 12)), col="lightsteelblue")
  
  #tetraploidy vs trisomy etc in suspected non-tetraploidy sets
  plot(main="tumours with > 10% diploid", nonTetraCNVs$tetraploid, nonTetraCNVs$triploid, xlim = c(0,1), ylim=c(0,1), pch=20, col="red", xlab = "proportion of genome at tetraploidy", ylab = "proportion of genome")
  points(nonTetraCNVs$tetraploid, nonTetraCNVs$Bloss, col="blue", pch=20)
  points(nonTetraCNVs$tetraploid, nonTetraCNVs$loss, col="lightsteelblue", pch=20)
  
  #fit regression lines
  R2fit1 <- lm(triploid ~ tetraploid, data=nonTetraCNVs)
  abline(R2fit1, col="red")
  R2fit2 <- lm(Bloss ~ tetraploid, data=nonTetraCNVs)
  abline(R2fit2, col="blue")
  R2fit3 <- lm(loss ~ tetraploid, data=nonTetraCNVs)
  abline(R2fit3, col="lightsteelblue")
  
  text(x=0.6, y=0.9, paste("R2:", format(summary(R2fit1)$r.squared, digits=4), "p=", round(summary(R2fit1)$coefficients[2,4], 12)), col="red")
  text(x=0.6, y=0.8, paste("R2:", format(summary(R2fit2)$r.squared, digits=4), "p=", round(summary(R2fit2)$coefficients[2,4], 12)), col="blue")
  text(x=0.6, y=0.7, paste("R2:", format(summary(R2fit3)$r.squared, digits=4), "p=", round(summary(R2fit3)$coefficients[2,4], 12)), col="lightsteelblue")
  
dev.off()




########### other plotting functions ########### 

#plot ploidy percentages (total samples)
pdf(file=paste(sampleList[1,6], plotDir, "ploidyPlot.pdf", sep=""), onefile=TRUE, width=8, height=10)
par(xpd=TRUE, mar=c(4,6,2,2))

#start plotting
plot(1, 1, col="white", axes=F, xlim=c(0,5.5), ylim=c(0, nrow(sampleList)), xlab="", ylab="", main="")

#colour list
colList <- c("grey", "red", "darkorange", "goldenrod", "orchid")

#for each sample set plot rectangles
for(currPlot in 1:nrow(CNVTable)){
  if(CNVTable[currPlot, 2] == "mean"){
    colTemp <- paste(colList, "3", sep="")
  }else{
    colTemp <- colList
  }
  
  #dip rect
  rect(xleft = 0, xright = CNVTable[currPlot, 3], ybottom = currPlot, ytop = currPlot+1, col=colTemp[1], border = "white")
  
  #poly rect
  rect(xleft = 1.1, xright = 1.1+CNVTable[currPlot, 4], ybottom = currPlot, ytop = currPlot+1, col=colTemp[2], border = "white")
  
  #tetra rect
  rect(xleft = 2.2, xright = 2.2+CNVTable[currPlot, 5], ybottom = currPlot, ytop = currPlot+1, col=colTemp[3], border = "white")
  
  #tri rect
  rect(xleft = 3.3, xright = 3.3+CNVTable[currPlot, 6], ybottom = currPlot, ytop = currPlot+1, col=colTemp[4], border = "white")
  
  #Bloh rect
  rect(xleft = 4.4, xright = 4.4+CNVTable[currPlot, 7], ybottom = currPlot, ytop = currPlot+1, col=colTemp[5], border = "white")
  
  #add average ploidy value
  if(!is.na(CNVTable[currPlot, "av.ploidy"])){
    text(x=5.5, y=(currPlot+0.5), label=round(CNVTable[currPlot, "av.ploidy"], 3), pos = 4, cex = 0.5)
  
  }
  
  if(CNVTable[currPlot, 2] == "mean"){
    text(x=round(CNVTable[currPlot, 3], 3), y=(currPlot+0.5), label=round(CNVTable[currPlot, 3], 3), pos = 4, cex = 0.5)
    text(x=(1.1 + round(CNVTable[currPlot, 4], 3)), y=(currPlot+0.5), label=round(CNVTable[currPlot, 4], 3), pos = 4, cex = 0.5)
    text(x=(2.2 + round(CNVTable[currPlot, 5], 3)), y=(currPlot+0.5), label=round(CNVTable[currPlot, 5], 3), pos = 4, cex = 0.5)
    text(x=(3.3 + round(CNVTable[currPlot, 6], 3)), y=(currPlot+0.5), label=round(CNVTable[currPlot, 6], 3), pos = 4, cex = 0.5)
    text(x=(4.4 + round(CNVTable[currPlot, 7], 3)), y=(currPlot+0.5), label=round(CNVTable[currPlot, 7], 3), pos = 4, cex = 0.5)
  }
}  

#add proportion axis
axis(side=1, at=seq(0,1,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.6, las=2)
axis(side=1, at=seq(1.1,2.1,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.6, las=2)
axis(side=1, at=seq(2.2,3.2,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.6, las=2)
axis(side=1, at=seq(3.3,4.3,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.6, las=2)
axis(side=1, at=seq(4.4,5.4,0.25), labels=c("", "0.25", "0.5", "0.75", ""), cex.axis=0.6, las=2)

#add names
axis(side=2, at=seq(1.5,(nrow(CNVTable)+1), 1), labels=CNVTable[[2]], cex.axis=0.4, las=2, tick = FALSE, line=-1.8)

#add ploidy names
axis(side=3, at=seq(0.5, 5.5, 1), labels=c("diploid", "polyploid", "tetraploid", "triploid", "LOH", "av.ploidy"), cex.axis=0.8, tick = FALSE, line=-1.5)

#add lines
lines(x=c(1.1,1.1),  y=c(0, (nrow(CNVTable)+1)), col="black")
lines(x=c(2.2,2.2),  y=c(0, (nrow(CNVTable)+1)), col="black")
lines(x=c(3.3,3.3),  y=c(0, (nrow(CNVTable)+1)), col="black")
lines(x=c(4.4,4.4),  y=c(0, (nrow(CNVTable)+1)), col="black")
lines(x=c(0,0),  y=c(0, (nrow(CNVTable)+1)), col="black")

dev.off()