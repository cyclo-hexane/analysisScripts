# script plots graphs of kmer sample pick against sample number and proportion gained
# 1.uses sample sampleList as script input
# 2.determines proportion and raw counts of variants called from conformance total by kmer picks.


sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[1])
sampleNames <- subset(sampleNames, sampleNames[1]!="Polyp.08")

#load driver gene list
driverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table5-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)

standardCols <- c("region", "gene","CHROM", "POS", "pos2", "REF", "ALT")

totalMutat <- c()
mainList <- as.list(0)
for(x in 1:nrow(sampleNames)){
	subList <- subset(sampleList, V1==sampleNames[x,1])
	
	#input sample vct (.txt) file
	dataSet <- read.table(file=paste(subList[1,6], "1.platypusCalls/somaticTotal.0.05/", subList[1,1], "/", subList[1,1], ".snv.annoVar.variant_function.0.05.txt", sep=""), sep="\t", stringsAsFactors=FALSE, header=FALSE)
	
  #filter for non-synonymous variants
	#dataSet <- dataSet[dataSet[2]=="nonsynonymous SNV" | dataSet[2]=="stopgain", ]
  
  #remove forst row (hash out for genome data)
	dataSet <- dataSet[-1]
	
	noSamples <- subList[1,8]
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
	noSamples <- 4
  
  #if samples are over 4, sub-sample
  if(length(namesList) > 4){
    subNames <- sample(namesList, 4)
    dataSet <- dataSet[c(standardCols, paste(subNames, ".NV", sep=""), subNames)]
  }else{
    subNames <- namesList
  }
  #remove variants with zeros representation in these four samples
  remList <- c(0)
  counter <- 1
  for(currRow in 1:nrow(dataSet)){
    if(sum(dataSet[currRow, c((8+noSamples):ncol(dataSet))]) == 0){
      remList[counter] <- currRow
      counter <- counter + 1
    }
  }
  if(remList[1]!=0){
    dataSet <- dataSet[-remList,]
  }
  
  #make every pairwise combination of samples and store in list
	combinations <- as.list(0)
	for(y in 1:4){
		combinations[[y]] <- combn(subNames, y)
	}
	
	#make and populate results table with variants from dataSet 
	variantNumbers <- as.list(0)
	for(z in 1:4){
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
  
	counter<-0
	asStrings <- as.list(0)
	for(k in seq(1,(2*length(variantNumbers)),2)){
		counter <- counter+1
		asStrings[[k]] <- paste("samples taken: ", counter, sep="")
		asStrings[[k+1]] <- variantNumbers[[counter]][[1]]
	}
	lapply(asStrings, write, paste(subList[1,6],"4.biopsyAnalysis/", subList[1,1], ".sampleRep.txt", sep=""), append=TRUE)
}

#################### now plot grid distributions ####################
resTab <- data.frame(matrix(NA, nrow=16, ncol=3))
names(resTab) <- c("type", "noSamples", "minCapture")
pdf(file=paste(subList[1,6],"4.biopsyAnalysis/allSets.sampleRep.pdf", sep=""), width=10, height=10)
	par(mfrow=c(4,4), xpd=TRUE, mar=c(4,4,4,4))
	for(f in 1:nrow(sampleNames)){
		subList <- subset(sampleList, V1==sampleNames[f,1])
		subMainList <- mainList[[f]]
			
		noSamples <- 4
	
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
		minLine <- median(subPlotTab[!is.na(subPlotTab[[1]]), 1])
		percLine <- minLine / subPlotTab[!is.na(subPlotTab[[4]]), 4]
		
		minLine2 <- median(subPlotTab[!is.na(subPlotTab[[2]]), 2])
		percLine2 <- minLine2 / subPlotTab[!is.na(subPlotTab[[4]]), 4]
		
		minLine3 <- median(subPlotTab[!is.na(subPlotTab[[3]]), 3])
		percLine3 <- minLine3 / subPlotTab[!is.na(subPlotTab[[4]]), 4]

		#finally plot data
		boxplot(subPlotTab, col=c("red", "steelblue", "green"), main=subList[1,1], ylim=c(0, subPlotTab[1,4]), xlab="number of biopsies", ylab="number of variants", cex.main=1, names=c(1:(noSamples)), cex.axis=1, cex.lab=1)
		axis(side=4, at=seq(0, subPlotTab[1,4], length.out=5), labels=seq(0,1,0.25))	
		
		#add one biopsy line
		lines(x=c(1, 1), y=c(0, minLine), col="red")
		lines(x=c(1, 5), y=c(minLine, minLine), col="red")
		text(x=1.5, y=minLine*0.8, labels = round(percLine, digits = 2), col="red")
		
		#add two biopsy line
		lines(x=c(2, 2), y=c(0, minLine2), col="steelblue")
		lines(x=c(2, 5), y=c(minLine2, minLine2), col="steelblue")
		text(x=2.5, y=minLine*0.8, labels = round(percLine2, digits = 2), col="steelblue")
		
		#add three biopsy line
		lines(x=c(3, 3), y=c(0, minLine3), col="green")
		lines(x=c(3, 5), y=c(minLine3, minLine3), col="green")
		text(x=3.5, y=minLine*0.8, labels = round(percLine3, digits = 2), col="green")
		
		resTab[f,2] <- noSamples
		resTab[f,3] <- round((sampleMin/sampleMax), 2)
	}
dev.off()


