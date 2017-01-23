#script inputs sample list and driver list and calculates frequencies
makeSig <- function(tableIn){
  sigTable <- table(tableIn[4])
  #print(sampleList[z,1])
  #print(sigTable)
  mutSig <- matrix(0, nrow=12, ncol=1)
  rownames(mutSig) <- c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G")
  for(y in 1:12){
    if(rownames(mutSig)[y] %in% names(sigTable)){
      mutSig[y,1] <- sigTable[names(sigTable)==rownames(mutSig)[y]]
    }else{mutSig[y,1] <- 0}
    
  }
  
  #merge comparative substitution columns
  mutSigFinal <- matrix(0, nrow=6, ncol=1)
  rownames(mutSigFinal) <- c( "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  mutSigFinal["C>T",] <- mutSig["C>T",] + mutSig["G>A",]
  mutSigFinal["C>G",] <- mutSig["C>G",] + mutSig["G>C",]
  mutSigFinal["C>A",] <- mutSig["C>A",] + mutSig["G>T",]
  mutSigFinal["T>A",] <- mutSig["A>T",] + mutSig["T>A",]
  mutSigFinal["T>C",] <- mutSig["A>G",] + mutSig["T>C",]
  mutSigFinal["T>G",] <- mutSig["A>C",] + mutSig["T>G",]
  
  return(mutSigFinal)
}



outDir <- "9.driverEvents/"
#outFileName <- "carcinoma"
#outFileName <- "lynch"
outFileName <- "adenoma"

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[3]])
#sampleNames <- sampleNames[c(1:11)]  
#sampleNames <- sampleNames[c(24:27)]
sampleNames <- sampleNames[c(12:16, 20:23)]

#get driver list
driverList <- read.table(file="~/PhD/CRCproject/10.finalFigures/supp.table.01-driverMutations.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
keepList <- c()
keepCounter <- 1
for(currRow in 1:nrow(driverList)){
  if(driverList[currRow, "Set"] %in% sampleNames){
    keepList[keepCounter] <- currRow
    keepCounter <- keepCounter + 1
  }
}

driverList <- driverList[keepList, ]

driverList <- driverList[driverList[["type"]]=="stopgain" | driverList[["type"]]=="nonsynonymous SNV", ]
driverList[4] <- paste(driverList[["ref"]], ">", driverList[["alt"]], sep="")

plotTables <- makeSig(driverList)

#plot histogram of driver frequency
pdf(file = paste(sampleList[1,6], outDir, outFileName, ".driverSigs.pdf", sep=""), width = 5, height = 3)
  barplot(height = plotTables, beside = TRUE, col=c("blue", "black", "red", "grey90", "green", "pink") ,las=2, main="Mutational signatures across sample branches and root", ylab="number of variants")
  legend('topright',legend=row.names(plotTables), col=c("blue", "black", "red", "grey90", "green", "pink"), lty=1, xpd=TRUE, cex = 0.5, lwd = 3)
dev.off()
  

#get indel types and plot


