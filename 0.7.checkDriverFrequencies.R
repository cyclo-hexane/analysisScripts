#script inputs sample list and driver list and calculates frequencies

outDir <- "9.driverEvents/"
outFileName <- "carcinoma.DriverFreq"
#outFileName <- "lynch.DriverFreq"
#outFileName <- "adenoma.DriverFreq"

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[3]])
sampleNames <- sampleNames[c(1:3,5:20)]  
#sampleNames <- sampleNames[c(24:27)]
#sampleNames <- sampleNames[c(12:16, 20:23)]

#get driver list
driverList <- read.table(file="Dropbox/MSeq-CRC/mseq.writeup/suppTables.final/table.S4-driverMutations.v2.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
for(currDriver in 1:nrow(driverList)){
  driverList[currDriver, 3] <- strsplit(driverList[currDriver, 3], split = ":")[[1]][1]
}

#get SNVs
#driverSNVList <- driverList[driverList[["type"]]=="stopgain" | driverList[["type"]]=="nonsynonymous SNV", ]

#get indels
#driverIndelList <- driverList[driverList[["type"]]!="stopgain" & driverList[["type"]]!="nonsynonymous SNV", ]

#populate driver table with stats
driverTab <- read.csv(file="/Users/cross01/Dropbox/MSeq-CRC/mseq.writeup/suppTables.final/table.S3-driverGeneList.csv", header = TRUE, stringsAsFactors = FALSE)
driverTab[2:20] <- NA
names(driverTab) <- c("gene", sampleNames)
for(currSam in 1:length(sampleNames)){
  subDriver <- driverList[driverList[["set"]]==sampleNames[currSam], ]
  for(currSNV in 1:nrow(subDriver)){
    driverTab[driverTab[["gene"]]==subDriver[currSNV, "gene"], sampleNames[currSam]] <- 1
  }
}

statsDriverTab <- data.frame(matrix(0, nrow = nrow(driverTab), ncol = 5))
names(statsDriverTab) <- c("gene", "freqAd", "freqCa", "pvalue",) 
statsDriverTab[1] <- driverTab[[1]]
adenomaDriverTab <- driverTab[c(1,12:20)]
carDriverTab <- driverTab[1:11]
for(currRow in 1:nrow(statsDriverTab)){
  tempGene <- statsDriverTab[currRow, 1]
  statsDriverTab[currRow, "freqAd"] <- sum(adenomaDriverTab[adenomaDriverTab[[1]]==tempGene, c(2:10)], na.rm = TRUE)
  statsDriverTab[currRow, "freqCa"] <- sum(carDriverTab[carDriverTab[[1]]==tempGene, c(2:10)], na.rm = TRUE)
  
  if(statsDriverTab[currRow, "freqCa"] > 1 | statsDriverTab[currRow, "freqAd"] > 1){
    statsTemp <- fisher.test(x = matrix(c(statsDriverTab[currRow, "freqAd"], statsDriverTab[currRow, "freqCa"], 9, 10), nrow = 2, ncol = 2, byrow = TRUE))
    statsDriverTab[currRow, "pvalue"] <- statsTemp$p.value
  }
}

write.csv(x = statsDriverTab, file="Dropbox/MSeq-CRC/mseq.writeup/suppTables.final/table.S5-driverEnrichment.v2.csv", quote = FALSE)
driverSNVList <- driverList

#make summary table
varTab <- data.frame(matrix(0, nrow=length(unique(driverList[[2]])), ncol=length(sampleNames)))
names(varTab) <- sampleNames
row.names(varTab) <- unique(driverList[[2]])

#loop through each sample and populate 
for(currSam in 1:length(sampleNames)){
  driverSub <- driverList[driverList[["Set"]]==sampleNames[currSam], ]
  
  currGenes <- unique(driverSub[["gene"]])
  
  for(add in 1:length(currGenes)){
    if("HIGH" %in% driverSub[driverSub[["gene"]]==currGenes[add], "cat"]){
      varTab[rownames(varTab)==currGenes[add], currSam] <- "H"
    }else{
      varTab[rownames(varTab)==currGenes[add], currSam] <- "L"
    }
  }
  
}

#add row totals
varTab[ncol(varTab)+1] <- 0
names(varTab)[ncol(varTab)] <- "counts"
varTab[ncol(varTab)+1] <- 0
names(varTab)[ncol(varTab)] <- "perc"
varTab[ncol(varTab)+1] <- 0
names(varTab)[ncol(varTab)] <- "Hcounts"
varTab[ncol(varTab)+1] <- 0
names(varTab)[ncol(varTab)] <- "Hperc"
for(currRow in 1:nrow(varTab)){
  tempTab <- table(varTab[currRow, sampleNames]=="L")
  if("TRUE" %in% names(tempTab)){
    varTab[currRow, "counts"] <- as.numeric(tempTab["TRUE"])
  }
  if(varTab[currRow, "counts"] != 0){
    varTab[currRow, "perc"] <- varTab[currRow, "counts"] / length(sampleNames)
  }
  
  tempTab <- table(varTab[currRow, sampleNames]=="H")
  if("TRUE" %in% names(tempTab)){
    varTab[currRow, "Hcounts"] <- as.numeric(tempTab["TRUE"])
  }
  if(varTab[currRow, "Hcounts"] != 0){
    varTab[currRow, "Hperc"] <- varTab[currRow, "Hcounts"] / length(sampleNames)
  }
}

#total up rows
varTab["total"] <- varTab[["counts"]] + varTab[["Hcounts"]]

#sort
varTab <- varTab[order(rownames(varTab), decreasing = FALSE), ]
varTab <- varTab[order(varTab[[ncol(varTab)]], decreasing = TRUE), ]


#plot histogram of driver frequency
pdf(file = paste(sampleList[1,6], outDir, outFileName, ".pdf", sep=""), width = 5, height = 3)
  barplot(t(as.matrix(varTab[c("Hperc", "perc")])), col=c("red", "grey60"), beside = FALSE, names.arg = rownames(varTab), las=2, cex.axis = 0.7, ylim=c(0,1), cex.names = 0.5)
dev.off()

#output driver query list
fileOut <- paste(sampleList[1,6], outDir, outFileName, ".txt", sep="")
write.table(varTab, file=fileOut, sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)  


