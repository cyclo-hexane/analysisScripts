# plot summary of CNV disruptions from cloneHD segmentations

##########################   notes   ##########################
#
# 
#
########################## libraries  ##########################

######################### subroutines ##########################


######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

setNames <- unique(sampleList[[1]])
setNames <- setNames[-c(4,12:24)]

#working dir
CNVholdingDir <- "7.CNVcalls.final/CNVannotation/"
CNVfileNames <- ".CNV.anno.txt"

#modelled dir
modHoldingDir <- "8.CNVphylogenetics/CNVphylo/"
modFileNames <- ".phyloCNVs.csv"

#clonal table
clonalTabSum <- data.frame(matrix(NA, nrow = 3, ncol = 30))
row.names(clonalTabSum) <- c("clonal", "regional", "unique")
names(clonalTabSum) <- c(setNames, paste(setNames, ".mod", sep=""), paste(setNames, ".space", sep=""))
clonalTabSum <- clonalTabSum[order(names(clonalTabSum))]

#size table
sizeTabSum <- data.frame(matrix(NA, nrow = 200, ncol = 30))
names(sizeTabSum) <- c(setNames, paste(setNames, ".mod", sep=""), paste(setNames, ".space", sep=""))
sizeTabSum <- sizeTabSum[order(names(sizeTabSum))]
pValues <- c(1:10)
pValuesClone <- c(1:10)
names(pValues) <- setNames

averageSeg <- c()
averageModSeg <- c()

for(currSam in 1:length(setNames)){
  #subset main list
  subSample <- subset(sampleList, sampleList[1]==setNames[currSam])
  
  sampleNames <- subSample[[2]]
  sampleNames <- sampleNames[-(subSample[1,7]+1)]
  noSamples <- length(sampleNames)
  
  #get modelled file
  modCNV <- read.csv(file=paste(subSample[1,6], modHoldingDir, subSample[1,1], modFileNames, sep=""), header=TRUE, stringsAsFactors=FALSE)
  modCNV <- modCNV[modCNV[["model"]]!="cannotModel",]
  modChrom <- unique(modCNV[[2]])
  clonalTabMod <- table(modCNV[["phyloLoc"]])
  totalClonalMod <- sum(clonalTabMod)
  #clonalTabMod <- clonalTabMod / totalClonalMod
  
  #get seg file
  segCNV <- read.table(file=paste(subSample[1,6], CNVholdingDir, subSample[1,1], CNVfileNames, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
  segCNV <- segCNV[segCNV[["chr"]]%in%modChrom, ]
  
  segCNV[segCNV[["phyloLoc"]]=="BiT", "phyloLoc"] <- "B"
  segCNV[segCNV[["phyloLoc"]]=="BiB", "phyloLoc"] <- "B"
  segCNV[segCNV[["phyloLoc"]]=="unknown >2 bifurcations", "phyloLoc"] <- NA
  segCNV <- segCNV[!is.na(segCNV[["phyloLoc"]]), ]
  clonalTab <- table(segCNV[["phyloLoc"]])
  totalClonal <- sum(clonalTab)
  #clonalTab <- clonalTab / totalClonal
  
  
  
  clonalTabSum[, subSample[1,1]] <- clonalTab[c("T", "B", "L")]
  clonalTabSum[, paste(subSample[1,1], ".mod", sep="")] <- clonalTabMod[c("T", "B", "L")]
  clonalTabSum[, paste(subSample[1,1], ".space", sep="")] <- c(0,0,0)
  
  averageSeg[currSam] <- mean(segCNV[["nloci"]])
  averageModSeg[currSam] <- mean(modCNV[["nloci"]])
  
  sizeTabSum[1:nrow(segCNV), subSample[1,1]] <- segCNV[["nloci"]]
  sizeTabSum[1:nrow(modCNV), paste(subSample[1,1], ".mod", sep="")] <- modCNV[["nloci"]]
  pValues[currSam] <- wilcox.test(segCNV[["nloci"]], modCNV[["nloci"]])$p.value
  
  #if("T" %in% names(clonalTab)){
  #  pValuesClone[currSam] <- binom.test(x = round(clonalTabMod["T"]*100, digits = 0), n = 100, p = clonalTab["T"])$p.value
  #}else{
  #  pValuesClone[currSam] <- 0
  #}
  
}
clonalTabSum[is.na(clonalTabSum)] <- 0
  
#plot clonal proportions
pdf(file=paste(subSample[1,6], "8.CNVphylogenetics/modDiff.pdf", sep=""), width=5, height=5)
  par(xpd=TRUE)
  barplot(as.matrix(clonalTabSum), las=2)
  text(x = seq(1,40,length.out = 10), y = 1, labels = pValuesClone, cex = 0.5, srt=45)
dev.off()

#plot sizes as boxplots
pdf(file=paste(subSample[1,6], "8.CNVphylogenetics/CNAMPsizeDiff.pdf", sep=""), width=5, height=5)
  boxplot(as.matrix(sizeTabSum), las=2, pch=20)
  text(x = seq(1,25,length.out = 10), y = 200, labels = pValues, cex = 0.5, srt=45)
dev.off()


mean(averageSeg)
mean(averageModSeg)

