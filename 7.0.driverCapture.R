# script plots graphs of kmer sample pick against sample number and proportion gained
# 1.uses sample sampleList as script input
# 2.determines proportion and raw counts of variants called from conformance total by kmer picks.


sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)

sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[-15]

#load driver gene list
driverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table3-driverDetails.csv", header=TRUE, stringsAsFactors=FALSE)

standardCols <- c("region", "gene","CHROM", "POS", "pos2", "REF", "ALT")
colVector <- c(rep("purple", 11), rep("red", 5))

totalMutat <- c()
mainList <- as.list(0)
for(x in 1:length(sampleNames)){
	subList <- subset(sampleList, V1==sampleNames[x])
	
	#input sample vct (.txt) file
	fileName <- paste(subList[1,6], "1.platypusCalls/somaticTotal.0.05/", subList[1,1], "/", subList[1,1], ".snv.annoVar.variant_function.0.05.txt", sep="")
	dataSet <- read.table(file=fileName, sep="\t", stringsAsFactors=FALSE, header=FALSE)
	
  #filter for non-synonymous variants
	#dataSet <- dataSet[dataSet[2]=="nonsynonymous SNV" | dataSet[2]=="stopgain", ]
  
  #remove forst row (hash out for genome data)
	dataSet <- dataSet[-1]
	
	noSamples <- (subList[1,8]-1)
	normalIndex <- subList[1,7]+1
  namesListNor <- subList[,2]
	normalName <- namesListNor[normalIndex] 
	namesListNV <- paste(namesListNor, ".NV", sep="")
	
  #name columns and remove normal
	names(dataSet) <- c(standardCols, namesListNV, namesListNor)
	remNormal <- c(normalName, paste(normalName, ".NV", sep=""))
	dataSet <- dataSet[which(remNormal[1]!=names(dataSet))]
	dataSet <- dataSet[which(remNormal[2]!=names(dataSet))]
  namesList <- namesListNor[-(normalIndex)]
  
  #make every pairwise combination of samples and store in list
	combinations <- as.list(0)
	for(y in 1:noSamples){
		combinations[[y]] <- combn(namesList, y)
	}
	
	#make and populate results table with variants from dataSet 
	variantNumbers <- as.list(0)
	for(z in 1:noSamples){
		variantNumbers[[z]] <- as.data.frame(0)
		for(w in 1:ncol(combinations[[z]])){
			subData <- subset(dataSet, select=c(combinations[[z]][,w]))
			if(ncol(subData)==1){
				variantNumbers[[z]][w,1] <- (length(subData[subData>0]))
			}else{
				variantNumbers[[z]][w,1] <- (nrow(subData[apply(subData, 1, function(x) !all(x==0)),]))
			}
		}
	}
	
	#save variant numbers to list
	mainList[[x]] <- variantNumbers
	totalMutat[x] <- nrow(dataSet)
	#totalMutat[x] <- as.numeric(table(dataSet[, namesList] > 0)["TRUE"])
  
# 	counter<-0
# 	asStrings <- as.list(0)
# 	for(k in seq(1,(2*length(variantNumbers)),2)){
# 		counter <- counter+1
# 		asStrings[[k]] <- paste("samples taken: ", counter, sep="")
# 		asStrings[[k+1]] <- variantNumbers[[counter]][[1]]
# 	}
# 	lapply(asStrings, write, paste(subList[1,6],"sampleRepresentation/", subList[1,1], ".sampleRep.txt", sep=""), append=TRUE)

}

#################### now plot grid distributions ####################
resTab <- data.frame(matrix(NA, nrow=16, ncol=3))
names(resTab) <- c("set", "noSamples", "minCapture")
pdf(file=paste(subList[1,6],"4.biopsyAnalysis/allSets.sampleRep.pdf", sep=""), width=10, height=10)
	par(mfrow=c(4,4), xpd=TRUE, mar=c(4,4,4,4))
	for(f in 1:length(sampleNames)){
		subList <- subset(sampleList, V1==sampleNames[f])
		subMainList <- mainList[[f]]
			
		noSamples <- (subList[1,8]-1)
	
		#get max row number for plotting tab below aswell as min and max values
		maxRowSub <- c(0)
		for(g in 1:length(subMainList)){
			maxRowSub[g] <- nrow(subMainList[[g]])
		}
		maxRowSub <- max(maxRowSub)
		
		#global max and min values
		sampleMax <- as.integer(subMainList[[(noSamples)]])
		sampleMin <- as.integer(min(subMainList[[1]]))
	
		#convert data to table with NAs inserted for plot function
		subPlotTab <- as.data.frame(matrix(data=NA, ncol=(noSamples), nrow=maxRowSub))
		for(k in 1:(noSamples)){
			tempList <- subMainList[[k]]
			if(nrow(tempList) == maxRowSub){
				subPlotTab[,k] <- tempList[[1]]
			}else{ tempList[(nrow(tempList)+1): maxRowSub,] <- NA
				subPlotTab[,k] <- tempList[[1]]}
		}
		
		#finally plot data
		boxplot(subPlotTab, ylim=c(0, sampleMax), main=paste("Sampling analysis ",subList[1,1], sep=""), xlab="number sampled", ylab="number of variants", cex.main=1, names=c(1:(noSamples)), cex.axis=1, cex.lab=1)
		axis(side=4, at=seq(0, sampleMax, length.out=5), labels=seq(0, 1, 0.25))	
		minCap <- min(subPlotTab[!is.na(subPlotTab[,1]) ,1])
		lines(x=c(1, 1), y=c(0, minCap), col="red", lwd=0.5)
		lines(x=c(1, (ncol(subPlotTab)+1)), y=c(minCap, minCap), col="red", lwd=0.5)
		resTab[f,1] <- subList[1,1]
		resTab[f,2] <- noSamples
		resTab[f,3] <- round((sampleMin/sampleMax), 2)
		
	}
dev.off()
minTab <- paste(subList[1,6],"4.biopsyAnalysis/1sampleCapture.txt", sep="")
write.csv(resTab, file=minTab, quote = FALSE, row.names=FALSE)



#################### assess driver capture from one sample ####################
driverPlotTab <- data.frame(matrix(NA, ncol=length(sampleNames), nrow=13))
names(driverPlotTab) <- sampleNames

driverNoVect <- c()

for(currSamp in 1:length(sampleNames)){
  subList <- subset(sampleList, V1==sampleNames[currSamp])
  
  #input sample vct (.txt) file
  dataSet <- read.table(file=paste(subList[1,6], "1.platypusCalls/nonSyn.0.05/", subList[1,1], "/", subList[1,1], ".snv.annoVar.exonic_variant_function.0.05.txt", sep=""), sep="\t", stringsAsFactors=FALSE, header=FALSE)
  
  #filter for non-synonymous variants
  dataSet <- dataSet[dataSet[2]=="nonsynonymous SNV" | dataSet[2]=="stopgain", ]
  
  #remove first row (hash out for genome data)
  dataSet <- dataSet[-1]
  
  noSamples <- subList[1,8]
  normalIndex <- subList[1,7]
  namesListNor <- subList[,2]
  namesListNV <- paste(namesListNor, ".NV", sep="")
  totalMutat <- nrow(dataSet)
  
  names(dataSet) <- c("region", "gene","CHROM", "POS", "pos2", "REF", "ALT", namesListNV, namesListNor)
  
  #get driver gene color scheme for plot
  dataSet[ncol(dataSet)+1] <- 0
  names(dataSet)[ncol(dataSet)] <- "driverFlag"
  for(n in 1:nrow(dataSet)){
    geneTemp <- strsplit(dataSet[n,2], ":")
    dataSet[n,2] <- geneTemp[[1]][1]
    if(dataSet[n,2] %in% driverList[[1]]){
      dataSet[n, ncol(dataSet)] <- 1
    }
  }
  
  driverTables <- subset(dataSet, dataSet["driverFlag"]==1)
  if(nrow(driverTables)!=0){
    totalDrivers <- nrow(driverTables)
    driverNoVect[currSamp] <- totalDrivers
    
    #remove normal sample
    namesList <- namesListNor[(normalIndex+1)]
    driverTables <- driverTables[,names(driverTables)!=c(namesList, paste(namesList, ".NV", sep=""))]
    namesList <- namesListNor[-(normalIndex+1)] 
    
    for(currDri in 1:length(namesList)){
      assessTable <- table(driverTables[namesList[currDri]] > 0)
      if("TRUE" %in% names(assessTable)){
        driverPlotTab[currDri, currSamp] <- as.numeric(assessTable["TRUE"]) / totalDrivers
      }else{
        driverPlotTab[currDri, currSamp] <- 0
      }
    }
  }else{
    driverNoVect[currSamp] <- 0
  }
  
} 

#get mean driver capture
meanDrivers <- c()
for(currSam in 1:length(sampleNames)){
  meanDrivers[currSam] <- mean(driverPlotTab[!is.na(driverPlotTab[currSam]), currSam])
}
colVector <- c(rep("purple", 11), rep("red", 5))
meanDrivers <- order(meanDrivers)
colVector <- colVector[meanDrivers]
driverPlotTab <- driverPlotTab[meanDrivers]

pdf(file=paste(subList[1,6],"4.biopsyAnalysis/driverCap.pdf", sep=""), width=8, height=6)
  par(xpd=TRUE, mar=c(8,5,5,5))
  #finally plot data
  boxplot(driverPlotTab, xaxt='n', main="driver mutation capture from one sample", xlab="", ylab="proportion of obverved driver mutations", cex.lab=1.5, col=colVector)
  axis(1, at=c(1:16), labels=names(driverPlotTab), las=2)
dev.off()

write.table(driverPlotTab, file="~/PhD/CRCproject/4.biopsyAnalysis/1samDriverCapture.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
