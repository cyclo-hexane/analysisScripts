# statistically assesses trunk and non-trunk CNV size
# additionally assesses distrubution of loss, gain CNVs

##########################   notes   ##########################

########################## libraries  ##########################

######################### subroutines ##########################

orderVariants <- function(tempTableSub, noSamInSet, samNames){
	#make blank column for ordering value
	tempTableSub[ncol(tempTableSub)+1] <- NA
	tempTableSub[ncol(tempTableSub)+1] <- NA
	
	for(k in 1:nrow(tempTableSub)){
		assessRow <- tempTableSub[k,9:(8+ noSamInSet)] !=0
		assessTable <- table(assessRow)
		
		#catch bulk variants not in set 
		if(names(assessTable)[1]=="FALSE" & as.integer(assessTable["FALSE"])== noSamInSet){
			tempTableSub[k,(ncol(tempTableSub)-1)] <- NA
			next
		}
		#if trunkal mark as 3
		if(as.integer(assessTable["TRUE"])== noSamInSet){
			tempTableSub[k,(ncol(tempTableSub)-1)] <- 3
		}
		if(as.integer(assessTable["TRUE"])==1){
			tempTableSub[k,(ncol(tempTableSub)-1)] <- 1
		}
		if(as.integer(assessTable["TRUE"])!=1 & as.integer(assessTable["TRUE"])!= noSamInSet){
			tempTableSub[k,(ncol(tempTableSub)-1)] <- 2
		}	
	}
	
	#remove bulk marked NA variants
	tempTableSub <- tempTableSub[!is.na(tempTableSub[ncol(tempTableSub)-1]), ]
	
	#now mark row AF sums for secondary ordering
	for(k in 1:nrow(tempTableSub)){
		assessRow <- tempTableSub[k,9:(8+ noSamInSet)]				
		tempTableSub[k,ncol(tempTableSub)] <- sum(assessRow)	
	}
	
	#final sort (clonal catagory and AF sum only)
	sortTable <- tempTableSub[ order( -tempTableSub[(ncol(tempTableSub)-1)], -tempTableSub[ncol(tempTableSub)] ), ]

	#name samples
	names(sortTable) <- c("varNo", "type", "gene", "chrom", "pos", "homoplasy", "ref", "alt", paste(samNames, "_AF", sep=""), "phyloLoc", "sumAF") 

	sortTable["chrom"] <- as.character(sortTable[["chrom"]])

	return(sortTable)
}

plotCNVMap <- function(ordDataSub, maxY, outputLocSub, driverListSub, statsSub){	
	sampleNamesSub <- names(ordDataSub)
	pdf(file=outputLocSub, width=20, height=10)
	par(xpd=TRUE, mar=c(10,3,3,3))
	plot(1, 1, col="white", axes=F, xlim=c(0,(length(ordDataSub)*5)), ylim=c(0,maxY) ,xlab="", ylab="", main="", cex.axis=2)
	#sample counter
	xcount <- 1	
	normalCounter <- 1
	#for each cancer sample plot bars
	for(cancer in 1:length(sampleNamesSub)){
		#get current sample Mb data and info
		currentSample <- ordDataSub[[cancer]]
		#get driver gene color scheme for plot
		driverCols <- data.frame(matrix(NA, ncol=1, nrow=nrow(currentSample)))
		driverCols[1] <- "black"
		for(n in 1:nrow(currentSample)){
			psplit1st <- strsplit(currentSample[n,3], split=":")
			for(currSplit in length(psplit1st[[1]])){
				psplit <- strsplit(psplit1st[[1]][currSplit], split="p")
				qsplit <- strsplit(psplit1st[[1]][currSplit], split="q")
				if(length(psplit[[1]]) > length(qsplit[[1]])){
					psplit <- paste(psplit[[1]][1], "p", sep="")
					if(psplit %in% driverListSub[[1]]){
						driverCols[n,1] <- "red"
					}			
				}else{
					qsplit <- paste(qsplit[[1]][1], "q", sep="")
					if(qsplit %in% driverListSub[[1]]){
						driverCols[n,1] <- "red"
					}
				}	
			}
		}

		#reset counters and table
		ycount <- 0
		ycountSNV <- 0
		textCounter <- 1
		textCoOrd <- data.frame(matrix(NA, ncol=1, nrow=nrow(currentSample)))
		#create rectangles on plot for current sample (CNV and SNV)
		for(currCNV in 1:nrow(currentSample)){
			#get bar colour (trunk, branch, leaf)		
			if(currentSample[currCNV,"phyloLoc"]=="T"){
				borderCol <- "steelblue"
			}
			if(currentSample[currCNV,"phyloLoc"]=="B"){
				borderCol <- "goldenrod"
			}
			if(currentSample[currCNV,"phyloLoc"]=="L"){
				borderCol <- "salmon"
			}
			#get Mb of current variant
			currentMb <- currentSample[currCNV,"Mb"]
			#plot Mb box and label for current sample
			rect(xleft=xcount, xright=xcount+1, ybottom=ycount, ytop=(ycount+ currentMb), col=borderCol, border="black")
			ycount <- ycount + currentMb + 0.01
			#get text co-ordinates
			textCoOrd[textCounter,1] <- ycount
			textCoOrd[textCounter,2] <- currentMb
			textCoOrd[textCounter,3] <- currentSample[currCNV,"cytoBand"]
			textCounter <- textCounter + 1
		}		

		#label up gene names in loop (to get different colors)
		textCoOrd[4] <- driverCols
		#filter labels for only large regions
		#textCoOrd <- subset(textCoOrd, textCoOrd[2]>20 | textCoOrd[4]=="red")
		textCoOrd <- subset(textCoOrd, textCoOrd[2]>50)
		if(nrow(textCoOrd) > 0){
			for(p in 1:nrow(textCoOrd)){
				text(y=textCoOrd[p,1], x=xcount+0.7, labels=textCoOrd[p,3], col=textCoOrd[p,4], cex=0.6, pos=4, srt=0)
			}	
		}
		xcount <- xcount + 5
	}
	#add y-axis and labels
	axis(2, at=seq(0,round(maxY, digits=-2), 250), line=-2, las=2)
	axis(1, at=seq(1.5,(length(sampleNames)*4.78), length=16), labels=sampleNamesSub, las=2, lwd=0, line=-2)
	text(y=-400, x=6, labels=paste("Fishers test, driver trunkal or non-trunkal: p-value =" ,round(statsSub, digits=2)))
	dev.off()
}

plotFunction <- function(resTableSub, outLocSub){
	resTableSub <- resTableSub[c(1:3),]
	fileName <- outLocSub
	pdf(file=fileName, height=8, width=8)
		par(mar=c(8,5,10,5), xpd=TRUE)
		barplot(as.matrix(resTableSub), col=c("steelblue", "goldenrod", "salmon"), main="variant distributions by set", ylab="", cex.lab=1.5, xaxt='n', space=1)
		axis(side=1, labels=names(resTable), at=seq(1.5,(length(sampleNames)*2),2), las=2)
		legend((length(names(resTable))*2)+1, y=600, legend=c("trunk", "branch", "leaf"), col=c("steelblue", "goldenrod", "salmon"), pch=15)
	dev.off()
}



######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)


#working dir
holdingDir <- "1.platypusCalls/nonSyn.0.05/"
namePrepended <- ".snv.annoVar.exonic_variant_function.0.05.txt"
CNVholdingDir <- "8.CNVphylogenetics/CNVphylo/"
CNVfileNames <- ".CNVsLoads.csv"
plotDir <- "2.phylogenetics/"

sampleList <- sampleList[sampleList[1]!="Polyp.08.WGS", ]

sampleNames <- unique(sampleList[[1]])
noSets <- length(sampleNames)

##################### 1.mutational loads for whole exome ##################### 

exomeTable <- data.frame(matrix(NA, ncol=length(sampleNames), nrow=3))
row.names(exomeTable) <- c("trunkal", "branch", "leaf")
names(exomeTable) <- sampleNames
 
for(j in 1:length(sampleNames)){
	print(paste("#### making table for sample ", sampleNames[j], " ####",sep=""))
	
	#subset main list
	subSample <- subset(sampleList, sampleList[1]==sampleNames[j])
	
	setName <- unique(subSample[[1]])
	noSamples <- subSample[1,8]
	normalIndex <- 1+subSample[1,7]
	samNames <- subSample[[2]]
	
	
	#setup input/output names
	dataIn <- read.table(file=paste(subSample[1,6], holdingDir, setName,"/", setName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)	
	tempTable <- dataIn
	
	#calculate allele frequecies
	for(y in 1:noSamples){
		tempTable[(ncol(tempTable)+1)] <- tempTable[(8+ noSamples +y)]/tempTable[(8 + y)]
	} 
		
	#remove unwanted columns	
	tempTable <- tempTable[,-c(9:(8+noSamples*2))]
	
	#remove normal
	tempTable <- tempTable[,-(8+normalIndex)]
	noSamples <- noSamples-1
	samNames <- samNames[-normalIndex]
	
	sortedTable <- orderVariants(tempTable, noSamples, samNames)
	
	variantCounts <- table(sortedTable["phyloLoc"])
	tempSum <- 100 / sum(variantCounts)
  
	exomeTable[1,j] <- variantCounts["3"] * tempSum
	exomeTable[2,j] <- variantCounts["2"] * tempSum
	exomeTable[3,j] <- variantCounts["1"] * tempSum
	
}
fileOutName <- paste(subSample[1,6], holdingDir, "nonSynCounts.csv", sep="")
write.csv(exomeTable, file=fileOutName)


##################### 2.acquire CNV mutational counts ##########################


CNVTable <- data.frame(matrix(NA, ncol=length(sampleNames), nrow=3))
row.names(CNVTable) <- c("trunkal", "branch", "leaf")
names(CNVTable) <- sampleNames

for(j in 1:length(sampleNames)){
	print(paste("#### making CNV table for sample ", sampleNames[j], " ####",sep=""))
	
	#subset main list
	subSample <- subset(sampleList, sampleList[1]==sampleNames[j])
	
	setName <- unique(subSample[[1]])
	noSamples <- subSample[1,8]
	normalIndex <- 1+subSample[1,7]
	samNames <- subSample[[2]]
	
	#setup input/output names
	dataIn <- read.table(file=paste(subSample[1,6], CNVholdingDir, setName, "/", setName, CNVfileNames, sep=""), sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
	
	divTot <- 100 / sum(sum(dataIn["nTrunk"]), sum(dataIn["nBranch"]), sum(dataIn["nLeaf"]))
	CNVTable[1,j] <- sum(dataIn["nTrunk"]) * divTot
	CNVTable[2,j] <- sum(dataIn["nBranch"]) * divTot
	CNVTable[3,j] <- sum(dataIn["nLeaf"]) * divTot
}

#produce comparison plot and stats
names(CNVTable) <- paste(names(CNVTable), ".CNV", sep="")
mergedVars <- cbind(exomeTable, CNVTable)
mergedVars <- mergedVars[order(names(mergedVars))]

#perform stats for below (trunk vs non-trunk) ######## results seem incorrect
pvalueVect <- c()
colSeq <- as.list(NA)
colSeq[[1]] <- c(1:11)
colSeq[[2]] <- c(12:16)
statCounter <- 1
for(currStat in 1:2){
  for(currCalc in 1:3){
    #null = pairs are from a different distribution
    pvalueVect[statCounter] <- wilcox.test(as.numeric(exomeTable[currCalc, colSeq[[currStat]]]), as.numeric(CNVTable[currCalc, colSeq[[currStat]]]), paired = TRUE, exact = TRUE)$p.value
    statCounter <- statCounter + 1  
  }
}
statTab <- matrix(pvalueVect, byrow = TRUE, nrow=2, ncol=3)
colnames(statTab) <- c("trunk", "branch", "leaf")
rownames(statTab) <- c("cancer", "adenoma")

#plot proportion distributions
fileName <- paste(subSample[1,6], plotDir, "phyloProportionComparison.pdf", sep="")
colVect <- c("steelblue", "goldenrod", "salmon")
adeCan <- c("cancer", "adenoma")
trunBran <- c("trunk", "branch", "leaf")
pdf(file=fileName, height=4, width=10)
par(mfrow=c(1,3), mar=c(8,5,5,8), xpd=TRUE)
  plotCounter <- 1
  for(currStat in 1:2){
    for(currCalc in 1:3){
     titleString <- paste(adeCan[currStat], "variant proportions,", trunBran[currCalc])
     boxplot(as.numeric(exomeTable[currCalc, colSeq[[currStat]]]), as.numeric(CNVTable[currCalc, colSeq[[currStat]]]), main=titleString, cex.lab=1.5, cex.axis=1.5, xaxt="n", col = colVect[currCalc], ylab="proportion of mutations")
     axis(side=1, labels=c("SNV", "CNV"), at=c(1:2), las=2, lwd=0, line=0.5, cex.axis=2)
     plotCounter <- plotCounter + 1  
    }
  }
dev.off()

#write stats file
write.csv(statTab, file=paste(subSample[1,6], plotDir, "phyloProportionComparison.stats.csv", sep=""), )
