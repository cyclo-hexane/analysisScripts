# script calculates diversity from conformace .vcf (converted using annoVar to .txt)
# 1. gets samples to be analysied from list
# 2. uses subroutine to create divergency table
# 3. outputs table
# extra: plots mutational burden boxplot

# updated version

#### notes ####
# .txt (.vcf) Files listed in above have be in format:
# CHROM, POS, REF, ALT, GEN[x].NR, GEN[x].NV, 4 x Annotation fields.
# changes list is the same file used in program 1.processVCFsnpSift.R
# only used to get sample names (line 51) Example: ~/PhD/CRCproject/sampleLists/nameChangesList.total.txt

#libraries
library(beeswarm)
library(lsr)
library(mclust)

################### subroutines ####################

#subroutine to get divergency table
makeDiv <- function(div, noSamp, sampleNamesSub, totalVar){
	#make output table
	divMat <- matrix(data=NA, nrow=noSamp, ncol=noSamp)
	colnames(divMat) <- sampleNamesSub
	rownames(divMat) <- sampleNamesSub
	for(w in 1:ncol(divMat)){
		divMat[w,w] <- 0
	}
	
	#numerate combinations
	expGrid <- expand.grid(sampleNamesSub, sampleNamesSub)
	expGrid2 <- t(apply(expGrid, 1, sort))
	combin <- expGrid2[!duplicated(expGrid2),]
	combinations <- subset(combin, combin[,1]!=combin[,2])
	combinations <- as.data.frame(combinations)
	combinations[3] <- NA
	
	#tally and populate divMat
	for(y in 1:nrow(combinations)){
		zeroCounts <- 0
		divCounts <- 0
		for(x in 1:nrow(div)){
			if(div[x,as.character(combinations[y,1])] > 0){
				compar1 <- 1
			}else{
				compar1 <- 0
			}
			if(div[x,as.character(combinations[y,2])] > 0){
				compar2 <- 1
			}else{
				compar2 <- 0
			}
			if(compar1 == compar2){
				if(compar1 == 0){
					zeroCounts <- zeroCounts + 1	
				}
			}
			else{
				divCounts <- divCounts + 1
			}	
		}
		divMat[as.character(combinations[y,1]),as.character(combinations[y,2])] <- divCounts
		divMat[as.character(combinations[y,2]),as.character(combinations[y,1])] <- divCounts
		combinations[y,3] <- totalVar - zeroCounts
	}
	returnList <- as.list(0)
	returnList[[1]] <- divMat
	returnList[[2]] <- combinations
	#return divergency table
	return(returnList)
}


#subset out variants for X% in at aleast one sample
filterAF <- function(confDataSub, sampleLocsSub, xPer){
	subCommand <- as.character()
	for(comSub in sampleLocsSub){
		currentCom <- paste("confDataSub[",comSub,"] > ", xPer," | ", sep="")
		subCommand <- paste(subCommand, currentCom, collapse=" ")
	}
	#remove last character and filter using subCommand string
	subCommand <- substr(subCommand, 1, nchar(subCommand)-2)
	confDataSub <- confDataSub[eval(parse(text=subCommand)),]
	return(confDataSub)
}


#make plotting table routine (currently defunct)
makePlotTable <- function(plotSamples, divFiles){
	sampleNameList <- unique(plotSamples[[1]])
	
	#setup plotting table
	plotMat <- data.frame(matrix(data=NA, ncol=2, nrow=1))
	names(plotMat) <- c("sampleName", "divergence")
	
	#get each diversity set into a table
	for(j in 1:length(sampleNameList)){
		sampleListSub <- subset(plotSamples, plotSamples[1]==sampleNameList[j], select=V2)
		normalIndex <- 1+as.integer(unique(subset(plotSamples, plotSamples[1]==sampleNameList[j], select=V7)))
		
		#remove normal from sampleList
		sampleListSub <- as.data.frame(sampleListSub[-normalIndex,])
		
		#get all sample conbinations for current set
		expGrid <- expand.grid(sampleListSub[[1]], sampleListSub[[1]])
		expGrid2 <- t(apply(expGrid, 1, sort))
		combin <- expGrid2[!duplicated(expGrid2),]
		combinations <- subset(combin, combin[,1]!=combin[,2])
		
		#loop through divergency table and add values to plotting table
		tempFile <- read.table(file=divFiles[j], sep="\t", stringsAsFactors=FALSE)
		print(tempFile)
		tempTable <- data.frame(matrix(data=NA, ncol=2, nrow=1))
		names(tempTable) <- c("sampleName", "divergence")
		for(k in 1:nrow(combinations)){
			tempTable[1,2] <- tempFile[combinations[k,1], combinations[k,2]]
			tempTable[1,1] <- j
			plotMat <- rbind(plotMat, tempTable)
		}
			
	}
	#remove NA in first row
	plotMat <- plotMat[-1,]
	return(plotMat)
}



################### main program ####################

arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=1){
	stop("\usage: > Rscript diversityAnalysis.R < fileList.csv > < holding directory > < .prependedName.vcf > ")
}

#sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filt.csv", header=FALSE, stringsAsFactors=FALSE)

#holdingDir <- arguments[2]
holdingDir <- "platypusCalls/conformanceData/"

outputDir <- "diversityAnalysis/"

#vcfName <- arguments[3]
exomeName <- ".snv.annoVar.exonic_variant_function.txt"
vcfName <- ".annoVar.variant_function.txt"

#get samplenames list
sampleNames <- sampleList[c(1,2,7)]
sampleNamesList <- unique(sampleNames[1])

#remove unwanted columns and make unique set table
sampleList <- sampleList[-c(3:4)]
sampleList <- unique(sampleList[c(1,3:6)])

#list of binary tables for burden plot
listTabs <- as.list(NA)
totalSetNames <- unique(sampleList[[1]])

#excluse MSI+ caner
totalSetNames <- totalSetNames[-c(4)]
sampleList <- sampleList[-4,]

totalDivTab <- as.list(NA)
counter <- 1
#change loop to 1:nrow(sampleList) for all samples
for(x in 1:length(totalSetNames)){
	exomeFile <- paste(sampleList[x,3], holdingDir, sampleList[x,1],"/", sampleList[x,1], exomeName, sep="")
	#varFile <- paste(sampleList[x,3], holdingDir, sampleList[x,1],"/", sampleList[x,1], vcfName, sep="")
	
	#input .txt files
	#confData <- read.table(file=varFile, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
	confDataExome <- read.table(file=exomeFile, sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)

	#get sample names
	currentNames <- subset(sampleNames, sampleNames[1]==sampleList[x,1], select=V2)
	currentSet <- as.character(unique(subset(sampleNames, sampleNames[1]==sampleList[x,1], select=V1)))

	#subset datasets
	#confData <- subset(confData, confData[1]!="exonic")
	#confData[1] <- "intronic"
	confDataExome <- confDataExome[,-1]
	
	noSamples <- sampleList[x,5]
	
	normalIndex <- 1+sampleList[x,4]
		
	confData <- confDataExome
	
	#calculate allele frequencies for both sets
	for(y in 1:noSamples){
		confData[(ncol(confData)+1)] <- confData[(7+noSamples+y)]/confData[7+y]
	}
	
	#remove unwanted columns	
	confData <- confData[,-c(8:(7+(noSamples*2)))]
	confData <- confData[-6]
	names(confData) <- c("type", "mut", "chrom", "pos", "ref", "alt", currentNames[[1]])
		
	#filter by AF (X% in at last one sample)
	#confData <- filterAF(confData, sampleLocs, 0)
	
	#subset tables for diversity analysis
	nonSynData <- subset(confData, confData[1]=="nonsynonymous SNV" | confData[1]=="stopgain")
	synData <- subset(confData, confData[1]=="synonymous SNV")
	
	#make divergency table from subroutine
	divTableTotalList <- makeDiv(confData, noSamples, currentNames[[1]], nrow(confData))
	divTableNonsynList <- makeDiv(nonSynData, noSamples, currentNames[[1]], nrow(nonSynData))
	divTableSynList <- makeDiv(synData, noSamples, currentNames[[1]], nrow(synData))
				
	#remove normals from tables and correct variables
	divTableTotal <- divTableTotalList[[1]]
	divTableNonsyn <- divTableNonsynList[[1]]
	divTableSyn <- divTableSynList[[1]]
	
	divTableTotal <- divTableTotal[-normalIndex, -normalIndex]
	divTableNonsyn <- divTableNonsyn[-normalIndex, -normalIndex]
	divTableSyn <- divTableSyn[-normalIndex, -normalIndex]
	noSamples <- noSamples-1
	currentNames <- currentNames[-normalIndex,]
	divList <- as.list(NA)
	divList[[1]] <- divTableTotal
	divList[[2]] <- divTableNonsyn
	divList[[3]] <- divTableSyn
		
	#make plotting table for each diversity table
	plotTableTemp <- data.frame(matrix(data=NA, ncol=4, nrow=1500))
	names(plotTableTemp) <- c("total", "nonSyn", "synon", "variants")
	for(divT in c(1:3)){
		currentTable <- divList[[divT]]
		dataTemp <- as.list(NA)
		nameTemp <- as.list(NA)
		for(samCol in 1:noSamples){
			 dataTemp[[samCol]] <- currentTable[-samCol,samCol]
			 nameTemp[[samCol]] <- paste(colnames(currentTable)[samCol], names(currentTable[-samCol,samCol]), sep=":")
		}
		plotTableTemp[[divT]][1:length(unlist(dataTemp))] <- unlist(dataTemp)
	}
	plotTableTemp[[4]][1:length(unlist(nameTemp))] <- unlist(nameTemp)
	
	#remove blank rows and unwanted columns
	plotTableTemp <- plotTableTemp[complete.cases(plotTableTemp),]	
	plotTableTemp <- plotTableTemp[duplicated(plotTableTemp[1:3]), ]
	plotTableTemp[5] <- NA
	tempStrings <- strsplit(plotTableTemp[[4]], split=":")
	for(currSplit in 1:length(tempStrings)){
		plotTableTemp[currSplit,4] <- tempStrings[[currSplit]][1]
		plotTableTemp[currSplit,5] <- tempStrings[[currSplit]][2]
	}

	#get variant numbers	
	varCombTemp <- divTableTotalList[[2]]
	
	for(currRow in 1:nrow(plotTableTemp)){
		plotTableTemp[currRow, 4] <- varCombTemp[(plotTableTemp[currRow, 4] == varCombTemp[[1]] | plotTableTemp[currRow, 4] == varCombTemp[[2]]) & (plotTableTemp[currRow, 5] == varCombTemp[[1]] | plotTableTemp[currRow, 5] == varCombTemp[[2]]), 3]
	}
	plotTableTemp <- plotTableTemp[-5]
	
	totalDivTab[[counter]] <- plotTableTemp
	listTabs[[counter]] <- plotTableTemp
	counter <- counter+1
}


##################### now compile into plottable table (3 types below) #####################

#compile table for box plot (grouped by mutant type)
compTableTotal <- data.frame(matrix(NA, ncol=((length(totalDivTab) * 4)), nrow=1000))
compCounter <- 0
compLocList <- seq(1,(length(totalDivTab) * 4), length(totalDivTab))
namePlotList <- c("total", "nonSyn", "synon", "diff")
for(currComp in 1:length(totalDivTab)){
	for(changeComp in 1:4){
		NAlist <- rep(NA,1000 - nrow(totalDivTab[[currComp]])) 
		compTableTotal[(compLocList[changeComp]+compCounter)] <- c(totalDivTab[[currComp]][[changeComp]], NAlist)
		names(compTableTotal)[ (compLocList[changeComp]+compCounter) ] <- paste(totalSetNames[currComp], "_", namePlotList[changeComp], sep="")
	}
	compCounter <- compCounter + 1
}
#plotting colours (for this table)
colVector <- c(rep("salmon3",4), rep("lightsalmon",5), rep("goldenrod",4), rep("lightgoldenrod",5), rep("steelblue",4), rep("lightsteelblue",5), rep("grey80",4), rep("grey50",5) )



#compile table for box plot (grouped by sample and diff removed with NA space inbetween samples)
compTableTotal <- data.frame(matrix(NA, ncol=((length(totalDivTab) * 4)), nrow=1000))
compCounter <- 0
namePlotList <- c("total", "nonSyn", "synon", "diff")
compCounter <- 1
for(currComp in 1:length(totalDivTab)){
	for(changeComp in 1:4){
		NAlist <- rep(NA,1000 - nrow(totalDivTab[[currComp]])) 
		compTableTotal[compCounter] <- c(totalDivTab[[currComp]][[changeComp]], NAlist)
		names(compTableTotal)[compCounter] <- paste(totalSetNames[currComp], "_", namePlotList[changeComp], sep="")
		if(changeComp == 4){
			compTableTotal[compCounter] <- NA
			names(compTableTotal)[compCounter] <- " "
		}
		compCounter <- compCounter + 1
	}
	
}
#plotting colours (for this table)
colVector1 <- c(rep(c("salmon3", "goldenrod", "steelblue", "grey80"), 11 ))
colVector2 <- c(rep(c("lightsalmon3", "lightgoldenrod", "lightsteelblue", "grey80"), 5 ))
colVector <- c(colVector1, colVector2)



#compile table for box plot with variant numbers on the y axis
totalRows <- 0
for(currIndex in 1:length(totalDivTab)){
	totalRows <- totalRows + nrow(totalDivTab[[currIndex]])
}
typeList <- c("red", "darkorange", "darkorange2", "red2", "orange2", "orange", "orangered", "orangered3", "red3", "maroon", "blue", "blue2", "cadetblue", "cadetblue3", "cyan")
compTableTotal <- data.frame(matrix(NA, ncol=3, nrow=totalRows))
compCounter <- 1
names(compTableTotal) <- c("varNo", "totalDiv", "type")
for(currComp in 1:length(totalDivTab)){
	#capture previous row number
	compPrev <- compCounter
	#get new row number
	compCounter <- (compPrev-1) + length(totalDivTab[[currComp]][[1]])
	#save values to plot table
	compTableTotal[c(compPrev:compCounter), 2] <- totalDivTab[[currComp]][[1]]
	compTableTotal[c(compPrev:compCounter), 1] <- totalDivTab[[currComp]][[4]]
	compTableTotal[c(compPrev:compCounter), 3] <- typeList[currComp]
	compCounter <- compCounter + 1
}
#get total distribution values
#compTableCancers <- compTableTotal[compTableTotal[[3]]=="red" ,]
#compTablePolyps <- compTableTotal[compTableTotal[[3]]=="blue" ,]

#cluster data
#clusterResults <- Mclust(compTableTotal[[2]])
#clustClass<- clusterResults$classification
#clustClass[which(clustClass==1)] <- "red"
#clustClass[which(clustClass==2)] <- "blue"

#now plot data
pdf(file="~/PhD/CRCproject/diversityAnalysis/divergencyBox.pdf", height=6, width=10)
	layout(matrix(c(1,2), nrow=2), widths=c(2.2,1), heights=c(0.5,3))
	par(mar=c(0.05,5,1.5,4))
	boxplot(compTablePolyps[[2]], compTableCancers[[2]], horizontal=TRUE, col=c("blue", "red"), xaxt='n', frame=FALSE, outline=FALSE, range=2, axes=FALSE)
	par(mar=c(5,5,1,4))
	plot(compTableTotal[[2]], compTableTotal[[1]], col=compTableTotal[[3]], ylim=c(50,300), pch=20, xlab="differing variants", ylab="total variants", cex.lab=1.5, cex.axis=1.5)
	legend("topright", legend=sampleList[[1]], col=typeList, pch=20, cex=0.5)
dev.off()

#stats
statRes <- ks.test(compTableCancers[[2]], compTablePolyps[[2]])


##################### plot one of the above tables #####################

#loop through samples and plot diversities
pdf(file=paste(sampleList[1,3], holdingDir, "divergencyTotal.pdf", sep=""), onefile=TRUE, width=10, height=10)
par(xpd=FALSE, mar=c(8,4.5,1,4.5), mfrow=c(2,1))
	plot.new()
	legend('bottomright',legend=c("total diversity", "non-syn diversity", "synon diversity"), col=c("salmon3", "goldenrod", "steelblue"), lty=1, xpd=TRUE, cex = 1, lwd = 5)
	boxplot(as.matrix(compTableTotal), col=colVector, las=2, xlab=" ", ylab=" ", cex=0.5, cex.lab=0.5)
	#text(seq(0.5, ncol(compTableTotal), length.out=36), labels=names(compTableTotal), srt=90)
	abline(v=seq(4,ncol(compTableTotal),4))
dev.off()

#### now perform statistics ####

setNames <- sampleList[[1]]
#setNames <- setNames[-c(2:5,9:11)]


##### cohens D, make results table #####
statsTable <- data.frame(matrix(data=NA, ncol=4, nrow=length(listTabs)))
names(statsTable) <- c("set", "mean-nonSyn", "mean-synon", "cohens-d")

for(counter in 1:length(listTabs)){
	statsTable[counter,1] <- setNames[counter]
	statsTable[counter,2] <- mean(listTabs[[counter]][[2]])
	statsTable[counter,3] <- mean(listTabs[[counter]][[3]])
	statTemp <- cohensD(x=listTabs[[counter]][[3]], y=listTabs[[counter]][[2]], method="paired")
	statsTable[counter,4] <- statTemp[1]
}
#write table
cohensTab <- paste(sampleList[1,3], holdingDir, "cohens.txt", sep="")
write.table(statsTable, file=cohensTab, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




##### perform kruskal wallis, non-parametric ANOVA on nonSyn-syn differences #####

#make results table
statsTable <- data.frame(matrix(data=NA, ncol=length(listTabs), nrow=length(listTabs)))
names(statsTable) <- setNames
row.names(statsTable) <- setNames

#perform tests
for(counter in 1:length(listTabs)){
	for(counter2 in 1:length(listTabs)){
		statsTemp <- wilcox.test(listTabs[[counter]][[4]], listTabs[[counter2]][[4]])
		statsTable[counter2,counter] <- statsTemp$p.value	
	}
}
compTab <- paste(sampleList[1,3], holdingDir, "stats.txt", sep="")
write.table(statsTable, file= compTab, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




#### perfrom linear regression on differences of nonSyn syn
#make cancer - adenoma list
canAdeList <- c(rep("cancer",11), rep("adenoma", 5))


#compile linear model table
linearRegTab <- data.frame(matrix(NA, nrow=1, ncol=6))
names(linearRegTab) <- c("total", "nonSyn", "synon", "diff", "variants", "type")
for(i in 1:length(listTabs)){
	listTemp <- listTabs[[i]]
	listTemp[6] <- canAdeList[i]
	names(listTemp)[6] <- "type"
	linearRegTab <- rbind(linearRegTab, listTemp)
}
#remove first column
linearRegTab <- linearRegTab[-1,]

#produce linear models
diffType <- summary(lm(data=linearRegTab, formula= diff ~ type))


#plot boxplot of differences
boxplot(linearRegTab[linearRegTab[[3]]=="cancer",4], linearRegTab[linearRegTab[[3]]=="adenoma",4] )

