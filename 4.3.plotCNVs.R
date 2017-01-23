# uses CNV disruptions from cloneHD segmentations

##########################   notes   ##########################
#
# plot summary of CNVs         
#

########################## subroutines  ##########################
determCol <- function(cnvTab){
  if(cnvTab[1, 6] == 1 & cnvTab[1, 7] == 1){
    plotCol <- "white"
  }else if(cnvTab[1, 6] == 0 & cnvTab[1, 7] == 1){
    plotCol <- "steelblue"
  }else if(cnvTab[1, 6] == 0 & cnvTab[1, 7] == 2){
    plotCol <- "darkblue"
  }else if(cnvTab[1, 6] == 1 & cnvTab[1, 7] == 2 | cnvTab[1, 6] == 0 & cnvTab[1, 7] == 3){
    plotCol <- "red"
  }else if(cnvTab[1, 6] == 2 & cnvTab[1, 7] == 2 | cnvTab[1, 6] == 1 & cnvTab[1, 7] == 3 | cnvTab[1, 6] == 0 & cnvTab[1, 7] == 4){
    plotCol <- "brown"
  }else{
    plotCol <- "purple"
  }
}


######################### libraries ##########################

######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

setNames <- unique(sampleList[[1]])
setNames <- setNames[-c(21:24)]

#working dir
CNVholdingDir <- "7.CNVcalls.final/CNVannotation/"
CNVfileNames <- ".CNV.anno.txt"

#chromosome positions
chromPos <- read.table(file="~/PhD/ReferenceGenome/chromosomeSizes.txt", stringsAsFactors = FALSE, header=FALSE, sep="\t")
chromPos[4] <- NA
names(chromPos) <- c("chrom", "start", "end", "cum")
for(currRow in 1:nrow(chromPos)){
  chromPos[currRow, 4] <- sum(as.numeric(chromPos[1:currRow, "end"]))
  if(currRow != 1){
    chromPos[currRow, 2] <- chromPos[(currRow-1), "cum"] + 1 
  }
} 
totalChromPos <- sum(as.numeric(chromPos[["end"]]))

CNVlist <- as.list(NA)

for(currSeg in 1:length(setNames)){
  print(paste("#### getting segs for sample ", setNames[currSeg], " ####",sep=""))
  setName <- setNames[currSeg]

  #set names for all samples in current set
  currentSamples <- sampleList[sampleList[[1]]==setName, ]
  samNames <- currentSamples[[2]]
  noSamples <- length(samNames)-1
  noNorNames <- samNames[-(currentSamples[1,7]+1)]
  
  #setup input/output names
  segCNVfile <- paste(sampleList[currSeg, 6], CNVholdingDir, setName, CNVfileNames, sep="")
  segCNV <- read.table(file=segCNVfile, sep="\t", fill=TRUE, stringsAsFactors=FALSE, header=TRUE)
  
  #segCNV["chr"] <- as.numeric(unlist(strsplit(segCNV[["chr"]], split = "chr"))[seq(2,(nrow(segCNV)*2),2)])
  #segCNV <- segCNV[order(segCNV[["chr"]]), ]

  CNVlist[[currSeg]] <- segCNV
  names(CNVlist)[[currSeg]] <- setName
  
}  




outputLoc <- paste(sampleList[currSeg, 6], CNVholdingDir, "CNV.plot.pdf", sep="")
pdf(file=outputLoc, width=8, height=5)
par(mar=c(2,2,2,2), xpd=FALSE, cex.main = 2.5) 
plot(1, 1, col="white", axes=F, xlim=c(0, totalChromPos), ylim=c(0, 101), xlab="chromosome positions", ylab="", main="")

currPlotRow <- 0

#plot segmentations for each sample
for(currSet in length(CNVlist):1){
  setName <- setNames[currSet]
  segCNV <- CNVlist[[currSet]]
  
  #set names for all samples in current set
  currentSamples <- sampleList[sampleList[[1]]==setName, ]
  samNames <- currentSamples[[2]]
  noSamples <- length(samNames)-1
  noNorNames <- samNames[-(currentSamples[1,7]+1)]
  
  for(currSam in 1:length(noNorNames)){
    currName <- noNorNames[currSam]
    if(setName == "Ha.22" | setName == "Ha.12" | setName == "Ha.6" | setName == "Ha.4"){
      currName <- paste("X", currName, sep="")
      
    }
    
    col1 <- paste(currName, "_Minor", sep="")
    col2 <- paste(currName, "_Major", sep="")
    
    currSegs <- segCNV[c("phyloLoc", "chr", "first.locus", "nloci", "last.locus", col1, col2)]
    
    #plot each seg
    for(currpl in 1:nrow(currSegs)){
      currSegInfo <- currSegs[currpl, ]
      
      #determine where to plot (which chromosome)
      chromTemp <- currSegInfo[1, "chr"]
      tempAdd <- chromPos[chromTemp, "start"]
      startPosTemp <- currSegInfo[1, "first.locus"] + tempAdd
      endPosTemp <- currSegInfo[1, "last.locus"] + tempAdd
      tempCol <- determCol(currSegInfo)
      rect(xleft = startPosTemp, xright = endPosTemp, ybottom = currPlotRow+0.05, ytop = currPlotRow+0.95, col = tempCol, border = NA)
    }
    #add biopsy name
    text(x = -10, y = currPlotRow+0.5, labels = currName, cex = 0.35)
    
    currPlotRow <- currPlotRow + 1
  }
  abline(h = currPlotRow)
}

#add chromosome lines
for(currLine in 1:24){
  abline(v = chromPos[currLine, "start"])
  text(x = (chromPos[currLine, "start"] + (chromPos[currLine, "end"] / 2)), y = currPlotRow+0.8, labels = currLine, cex = 0.5)
}

dev.off()




