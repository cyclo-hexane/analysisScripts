# make mutational load tables/plots for SNV and CNVs, both as counts
# and as a proportion (SNV and CNV plotted together for comparison)

######################## notes ########################

######################## libraries ########################


######################## subroutines ########################

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
  sampleNamesSub <- names(resTableSub)
  resTableSub <- resTableSub[c(1:3),]
	fileName <- outLocSub
	pdf(file=fileName, height=8, width=8)
		par(mar=c(8,5,10,5), xpd=TRUE)
		barplot(as.matrix(resTableSub), col=c("steelblue", "goldenrod", "salmon"), main="variant distributions by set", ylab="", cex.lab=1.5, xaxt='n', space=1)
		axis(side=1, labels=names(resTableSub), at=seq(1.5,(length(sampleNamesSub)*2),2), las=2)
		legend((length(names(resTableSub))*2)+1, y=600, legend=c("trunk", "branch", "leaf"), col=c("steelblue", "goldenrod", "salmon"), pch=15)
	dev.off()
}

######################## main program ########################


#input sampleList from commandline arguments
arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=4){
	stop("\n#### please use syntax > Rscript 3.2.mutationalLoads.R < sample list file > < holding directory > < prepended.name > < homoplasy tables > ####\n")
}
sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
holdingDir <- arguments[2]
namePrepended <- arguments[3]
holdingDir2 <- arguments[4]

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
holdingDir <- "1.platypusCalls/nonSyn.0.05/"
namePrepended <- ".snv.annoVar.exonic_variant_function.0.05.txt"
holdingDir3 <- "8.CNVphylogenetics/CNVphylo/"
CNVfileNames <- ".phyloCNVs.csv"
plotDir <- "10.finalFigures/"

#CNV driver list
driverListCNV <- read.csv(file="~/PhD/CRCproject/9.driverEvents/driverList.chromLoci.csv", header=FALSE, stringsAsFactors=FALSE)

#total sample list for filtering
#sampleListRem <- read.csv(file="~/PhD/CRCproject/masterSampleList.total.exclusions.csv", header=FALSE, stringsAsFactors=FALSE)

#cyto band table
cytoName <- "/Users/cross01/PhD/ReferenceGenome/cytoBand.txt"
cytoBands <- read.table(file=cytoName, sep="\t", header=FALSE, stringsAsFactors=FALSE)

sampleNames <- unique(sampleList[[1]])

#load driver gene list
driverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table5-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)



######################## assess SNV loads ########################
#outputs a resTable table of total variants

#setup variant results table
resTable <- data.frame(matrix(NA, ncol=length(sampleNames), nrow=3))
row.names(resTable) <- c("trunkal", "branch", "leaf")
names(resTable) <- sampleNames

#now process .vcf file and calculate mutational loads 
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
	
	#append row (for genome data only)
	#dataIn <- as.data.frame(append(dataIn, NA, after=1))
	
	#filter out synonymous variants
	#tempTable <- subset(dataIn , dataIn$V2=="nonsynonymous SNV" | dataIn$V2=="stopgain")
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
	
	resTable[1,j] <- variantCounts["3"]
	resTable[2,j] <- variantCounts["2"]
	resTable[3,j] <- variantCounts["1"]
	
}


######################## assess CNV loads ########################
#outputs a resCNVtable table of total Mb changed and plotCNVList of detailed loci intervals

#setup variant results table
resCNVtable <- data.frame(matrix(NA, ncol=length(sampleNames), nrow=3))
row.names(resCNVtable) <- c("trunkal", "branch", "leaf")
names(resCNVtable) <- sampleNames

plotCNVList <- as.list(NA)

maxYsub <- 0

#loop through each CNV file and get total Mb pairs disrupted
for(j in 1:length(sampleNames)){
	print(paste("#### making CNV table for sample ", sampleNames[j], " ####",sep=""))
	
	#subset main list
	subSample <- subset(sampleList, sampleList[1]==sampleNames[j])
	
	setName <- unique(subSample[[1]])
	noSamples <- subSample[1,8]-1
	normalIndex <- 1+subSample[1,7]
	samNames <- subSample[[2]]
	samNames <- samNames[-normalIndex]
	
	#setup input/output names
	dataIn <- read.table(file=paste(subSample[1,6], holdingDir3, setName, CNVfileNames, sep=""), sep=",", colClasses = c("character"), header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
	
	#input to main plotting table (for MB)
	totalMB <- sum(as.numeric(dataIn[["nloci"]]))
	resCNVtable["trunkal",j] <- sum(as.numeric(dataIn[dataIn[["phyloLoc"]]=="T", "nloci"])) / totalMB
	resCNVtable["branch",j] <- sum(as.numeric(dataIn[dataIn[["phyloLoc"]]=="B", "nloci"])) / totalMB
	resCNVtable["leaf",j] <- sum(as.numeric(dataIn[dataIn[["phyloLoc"]]=="L", "nloci"])) / totalMB
  
	#input to main plotting table (for counts)
	#resCNVtable["trunkal",j] <- nrow(dataIn[dataIn[["phyloLoc"]]=="T",])
	#resCNVtable["branch",j] <- nrow(dataIn[dataIn[["phyloLoc"]]=="B",])
	#resCNVtable["leaf",j] <- nrow(dataIn[dataIn[["phyloLoc"]]=="L",])
}


###################### reorder and normalize for plots ######################

#normalize SNV values (if needed)
resNorm <- resTable
for(currSamp in 1:length(sampleNames)){
	tempSum <- sum(resNorm[[currSamp]])
	for(currRow in 1:3){
		resNorm[currRow, currSamp] <- resNorm[currRow, currSamp] / tempSum
	}
}

#normalize CNV values (if needed)
resCNVnorm <- resCNVtable
for(currSamp in 1:length(sampleNames)){
	tempSum <- sum(resCNVnorm[[currSamp]])
	for(currRow in 1:3){
		resCNVnorm[currRow, currSamp] <- resCNVnorm[currRow, currSamp] / tempSum
	}
}

resCNVnorm <- resCNVnorm[order(resNorm[1,])]
resNorm <- resNorm[order(resNorm[1,])]
resCNVnorm[4,] <- resCNVnorm[2,] + resCNVnorm[3,]
resNorm[4,] <- resNorm[2,] + resNorm[3,]
resCNVnormStats <- round((resCNVnorm*100))
resNormStats <- round((resNorm*100))

#get total burden column and reorder
resTable[4,] <- 0
for(currRes in 1:length(sampleNames)){
	resTable[4, currRes] <- sum(resTable[, currRes])
}
row.names(resTable)[4] <- "total"
resCNVtable[4,] <- 0
for(currRes in 1:length(sampleNames)){
	resCNVtable[4, currRes] <- sum(resCNVtable[, currRes])
}
row.names(resCNVtable)[4] <- "total"
resCNVtable <- resCNVtable[order(resTable[4,])]
resTable <- resTable[order(resTable[4,])]


#merge CNV and SNV loads for correlation plot
mergedLoads <- data.frame(matrix(NA, nrow=length(sampleNames)*4, ncol=3))
names(mergedLoads) <- c("SNV", "CNV", "phylo")
rowCounter <- 1
for(currSam in 1:length(sampleNames)){
	for(currPhylo in 1:4){
		mergedLoads[rowCounter, 1] <- resTable[currPhylo, currSam]
		mergedLoads[rowCounter, 2] <- resCNVtable[currPhylo, currSam]
		mergedLoads[rowCounter, 3] <- row.names(resCNVtable[currPhylo,])
		rowCounter <- rowCounter + 1	
	}
}
#remove set.04
mergedLoads <- mergedLoads[-c(61:64),]

#get manually imput data order for plotCNVtable
reorderVect <- c(15, 12, 14, 16, 13, 4, 7, 2, 5, 10, 8, 1, 6, 9, 3, 11)

#reorder data CNV and SNV
plotCNVListCopy <- plotCNVList
namesVect <- names(resCNVtable)[reorderVect]
for(currReord in 1:length(namesVect)){
	plotCNVListCopy[[currReord]] <- plotCNVList[[namesVect[currReord]]]
}
names(plotCNVListCopy) <- namesVect


######################## ploting above ordered data ########################

#now SNV loads (re-ordering first) 
fileName <- paste(subSample[1,6], plotDir, "SNVloads.pdf", sep="")
resTable <- resTable[, c(3, 7, 9, 11, 12, 1, 2, 4, 5, 6, 8, 13, 14, 15, 16, 17, 10)]
resTable[17] <- NULL
plotFunction(resTable, fileName)


#plot CNV loads (re-ordering first) 
fileName <- paste(subSample[1,6], plotDir, "CNVloads.pdf", sep="")
#resCNVtable <- resCNVtable[, c(3, 7, 9, 11, 12, 1, 2, 4, 5, 6, 8, 13, 14, 15, 16, 17, 10)]
#resCNVtable[17] <- NULL
plotFunction(resCNVtable, fileName)

fileName <- paste(subSample[1,6], plotDir, "CNVclonalityTab.txt", sep="")
write.table(resCNVtable, file=fileName, sep="\t", quote = FALSE, row.names = TRUE)

#perform stats for below (trunk vs non-trunk)
statsList <- as.list(NA)
pvalueVect <- c(0)
for(currStat in 1:length(sampleNames)){
  statsTab <- data.frame(matrix(NA, nrow=2, ncol=2))
  names(statsTab) <- c("trunkal", "nonTrunkal")
  row.names(statsTab) <- c("SNV", "CNV")
  statsTab[1,1] <- resNormStats["trunkal",currStat]
  statsTab[1,2] <- resNormStats[4,currStat]
  statsTab[2,1] <- resCNVnormStats["trunkal",currStat]
  statsTab[2,2] <- resCNVnormStats[4,currStat]
  statsList[[currStat]] <- statsTab
  pvalueVect[currStat] <- chisq.test(statsTab)$p.value
}
names(pvalueVect) <- names(resNorm)

#plot merged proportions
resCNVnorm <- resCNVnorm[-c(2,3),]
resNorm <- resNorm[-c(2,3),]
names(resCNVnorm) <- paste(names(resCNVnorm), ".CNV", sep="")
mergedVars <- cbind(resNorm, resCNVnorm)
mergedVars <- mergedVars[order(names(mergedVars))]

#reorder p-values
pvalueVect <- pvalueVect[order(names(pvalueVect))]


fileName <- paste(subSample[1,6], plotDir, "totalVarplot.merged.pdf", sep="")
pdf(file=fileName, height=4, width=20)
par(mar=c(8,5,5,8), xpd=TRUE)
barplot(as.matrix(mergedVars), col=c("steelblue", "orange3"), main="normalized variant distributions by set", ylab="proportion of total burden", cex.lab=1.5, xaxt='n', space=c(2,0.5), ylim=c(0,1))
axis(side=1, labels=names(pvalueVect), at=seq(3, 75.5, length=length(sampleNames)), las=2, lwd=0, line=-0.5)
legend(x=75, y=1, legend=c("SNV trunk", "SNV non-trunk", "CNV trunk", "CNV non-trunk"), col=c("steelblue", "orange3", "lightsteelblue", "orange"), pch=15)

text(x=-3, y=1.05, labels=paste("chisq P-value:"))
for(currLab in 1:length(pvalueVect)){
  text(x=(currLab*4.5), y=1.05, labels=round(pvalueVect[currLab], digit=3)) 
}
dev.off()


#plot merged CNV and SNV loads for correlation analysis
mergedLoads[mergedLoads[[3]]=="trunkal", 3] <- "steelblue"
mergedLoads[mergedLoads[[3]]=="branch", 3] <- "goldenrod"
mergedLoads[mergedLoads[[3]]=="leaf", 3] <- "salmon"
mergedLoads[mergedLoads[[3]]=="total", 3] <- "grey40"
fileName <- paste(subSample[1,6], plotDir,"CNV-SNV.corr.pdf", sep="")
pdf(file=fileName, height=4, width=4)
	plot(mergedLoads[[2]], mergedLoads[[1]], xlab="total CNV load", ylab="total (N/S) SNV load", col=mergedLoads[[3]], pch=20)
dev.off()



#################### plot table of annotated/partitioned CNV and SNV values ####################


#get driver status for each CNV then perform stats for distribution
for(currSet in 1:length(plotCNVListCopy)){
	currentSample <- plotCNVListCopy[[currSet]]
	
	#get driver gene color scheme for plot
	driverCols <- data.frame(matrix(NA, ncol=1, nrow=nrow(currentSample)))
	driverCols[1] <- "benign"
	for(n in 1:nrow(currentSample)){
		psplit1st <- strsplit(currentSample[n,3], split=":")
		for(currSplit in length(psplit1st[[1]])){
			psplit <- strsplit(psplit1st[[1]][currSplit], split="p")
			qsplit <- strsplit(psplit1st[[1]][currSplit], split="q")
			if(length(psplit[[1]]) > length(qsplit[[1]])){
				psplit <- paste(psplit[[1]][1], "p", sep="")
				if(psplit %in% driverListCNV[[1]]){
					driverCols[n,1] <- "driver"
				}			
			}else{
				qsplit <- paste(qsplit[[1]][1], "q", sep="")
				if(qsplit %in% driverListCNV[[1]]){
					driverCols[n,1] <- "driver"
				}
			}	
		}
	}
	
	#add driver status to table	
	currentSample[5] <- driverCols[[1]]
	names(currentSample)[5] <- "status"
	plotCNVListCopy[[currSet]] <- currentSample	
}

#get stats table data
statsCNVtable <- data.frame(matrix(NA, nrow=2, ncol=2))
names(statsCNVtable) <- c("driver", "benign")
row.names(statsCNVtable) <- c("trunk", "nonTrunk")
trunkDriverTotals <- 0
trunkBenignTotals <- 0
nonTrunkDriverTotals <- 0
nonTrunkBenignTotals <- 0
for(currSet in 1:length(plotCNVListCopy)){
	currentSample <- plotCNVListCopy[[currSet]]
	trunkCNVs <- subset(currentSample, currentSample["phyloLoc"]=="T")
	nonTrunkCNVs <- subset(currentSample, currentSample["phyloLoc"]=="B" | currentSample["phyloLoc"]=="L")
	
	if(nrow(trunkCNVs) == 0){
		trunkDriverTotals <- 0
		trunkBenignTotals <- 0
	}else if(nrow(nonTrunkCNVs) == 0){
		nonTrunkDriverTotals <- 0
		nonTrunkBenignTotals <- 0
	}else{
		trunkDriver <- subset(trunkCNVs, trunkCNVs["status"]=="driver")
		trunkBenign <- subset(trunkCNVs, trunkCNVs["status"]=="benign")
		nonTrunkDriver <- subset(nonTrunkCNVs, nonTrunkCNVs["status"]=="driver")
		nonTrunkBenign <- subset(nonTrunkCNVs, nonTrunkCNVs["status"]=="benign")
		if(nrow(trunkDriver)!=0){
			trunkDriverTotals <- trunkDriverTotals + sum(trunkDriver["Mb"])
		}
		if(nrow(trunkBenign)!=0){
			trunkBenignTotals <- trunkBenignTotals + sum(trunkBenign["Mb"])
		}
		if(nrow(nonTrunkDriver)!=0){
			nonTrunkDriverTotals <- nonTrunkDriverTotals + sum(nonTrunkDriver["Mb"])
		}
		if(nrow(nonTrunkBenign)!=0){
			nonTrunkBenignTotals <- nonTrunkBenignTotals + sum(nonTrunkBenign["Mb"])	
		}
	}
}		
statsCNVtable[1,1] <- trunkDriverTotals
statsCNVtable[1,2] <- trunkBenignTotals
statsCNVtable[2,1] <- nonTrunkDriverTotals
statsCNVtable[2,2] <- nonTrunkBenignTotals

#normalize for calculations
trunkTot <- sum(statsCNVtable[1,])
nonTunkTot <- sum(statsCNVtable[2,])
statsCNVtable[1,1] <- statsCNVtable[1,1] / trunkTot
statsCNVtable[1,2] <- statsCNVtable[1,2] / trunkTot
statsCNVtable[2,1] <- statsCNVtable[2,1] / nonTunkTot
statsCNVtable[2,2] <- statsCNVtable[2,2] / nonTunkTot
statsCNVtable <- round(statsCNVtable*100, digits=0)
statsCNVResults <- fisher.test(statsCNVtable)


#plot CNV distributions
outputLoc <- paste(subSample[1,6], plotDir, "CNVloads.anno.pdf", sep="")
plotCNVMap(plotCNVListCopy, maxYsub, outputLoc, driverListCNV, statsCNVResults$p.value)
