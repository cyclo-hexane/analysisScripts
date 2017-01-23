# make mutational load tables/plots for SNV and CNVs, both as counts
# and as a proportion (SNV and CNV plotted together for comparison)

######################## notes ########################
# assesses SNV clonality by subsampling variants to 4 biopsies
# plots variants as boxplot with confidence interval
#
#
######################## libraries ########################


######################## subroutines ########################


######################## main program ########################


#input sampleList from commandline arguments
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- "1.platypusCalls/somaticTotal.0.01/"
plotDir <- holdingDir
namePrepended <- ".snv.annoVar.variant_function.0.01.txt"

sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[-c(4,20:24)]


######################## assess SNV clonality ########################

#plot table
plotTab <- data.frame(matrix(NA, ncol=length(sampleNames), nrow=3))
names(plotTab) <- sampleNames
row.names(plotTab) <- c("trunk", "branch", "leaf")

#clonal subsampling list
clonalList <- as.list(NA)

#calculate mutational proportions 
for(j in 1:length(sampleNames)){
	print(paste("#### making plot for sample ", sampleNames[j], " ####",sep=""))
	
	#subset main list
	subSample <- subset(sampleList, sampleList[1]==sampleNames[j])
	
	setName <- unique(subSample[[1]])
	noSamples <- subSample[1,8]
	normalIndex <- 1+subSample[1,7]
	samNames <- subSample[[2]]
	normalName <- samNames[normalIndex]
	samNamesNoNorm <- samNames[-normalIndex]
	
	#setup input/output names
	dataIn <- read.table(file=paste(subSample[1,6], holdingDir, setName,"/", setName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
	
	#name columns
  names(dataIn) <- c("line", "type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(samNames, ".NR", sep=""), samNames)
		
  #make clonality table
  clonalityTab <- data.frame(matrix(NA, ncol=3, nrow=100))
	names(clonalityTab) <- c("trunk", "branch", "leaf")
  
  #get subsample names list
	if(length(samNamesNoNorm) > 4){
	  
	  #get unique sample quartet-wise combinations
	  sampComb <- as.data.frame(combn(samNamesNoNorm, m = 4))
	  
	  #get clonality score for each combination
	  for(currComb in 1:ncol(sampComb)){
	    #current samples
	    samTemp <- as.character(sampComb[[currComb]])
	    
	    #get samples
	    tempDataIn <- dataIn[, samTemp]

	    #count variants
	    trunkCounter <- 0
	    branchCounter <- 0
	    leafCounter <- 0
	    for(currRow in 1:nrow(tempDataIn)){
	     varTemp <- table(tempDataIn[currRow, ] > 0)
	     if("TRUE" %in% names(varTemp)){
	       tempVarHolder <- as.numeric(varTemp["TRUE"])
	     }else{
	       #variant doesn't exist for this sample
	       next
	     }
	     
	     if(tempVarHolder == 4){
	       trunkCounter <- trunkCounter + 1
	     }else if(tempVarHolder == 1){
	       leafCounter <- leafCounter + 1
	     }else{
	       branchCounter <- branchCounter + 1
	     }
	     
	     #total variants
	     totalVar <- trunkCounter + branchCounter + leafCounter
	     
	     #add proportions to clonality table
	     clonalityTab[currComb, "trunk"] <- trunkCounter / totalVar
	     clonalityTab[currComb, "branch"] <- branchCounter / totalVar
	     clonalityTab[currComb, "leaf"] <- leafCounter / totalVar
	    }
	    
	    #store clonality table
	    clonalList[[j]] <- clonalityTab
	    
	    #add averages to plot tab
	    plotTab[1, j] <- median(clonalityTab[[1]], na.rm = TRUE)
	    plotTab[2, j] <- median(clonalityTab[[2]], na.rm = TRUE)
	    plotTab[3, j] <- median(clonalityTab[[3]], na.rm = TRUE)
	  }
	  
	}else{
	  #no need to subsample
	  
	  #count variants
	  trunkCounter <- 0
	  branchCounter <- 0
	  leafCounter <- 0
	  for(currRow in 1:nrow(dataIn)){
	    varTemp <- table(dataIn[currRow, samNamesNoNorm] > 0)
	    if("TRUE" %in% names(varTemp)){
	      tempVarHolder <- as.numeric(varTemp["TRUE"])
	    }else{
	      #variant doesn't exist for this sample
	      next
	    }
	    
	    if(tempVarHolder == length(samNamesNoNorm)){
	      trunkCounter <- trunkCounter + 1
	    }else if(tempVarHolder == 1){
	      leafCounter <- leafCounter + 1
	    }else{
	      branchCounter <- branchCounter + 1
	    }
	  } 
	  #total variants
	  totalVar <- trunkCounter + branchCounter + leafCounter
	  
	  #add proportions to clonality table
	  clonalityTab[1, "trunk"] <- trunkCounter / totalVar
	  clonalityTab[1, "branch"] <- branchCounter / totalVar
	  clonalityTab[1, "leaf"] <- leafCounter / totalVar
	  
	  #store clonality table
	  clonalList[[j]] <- clonalityTab
	  
	  #add averages to plot tab
	  plotTab[1, j] <- mean(clonalityTab[[1]], na.rm = TRUE)
	  plotTab[2, j] <- mean(clonalityTab[[2]], na.rm = TRUE)
	  plotTab[3, j] <- mean(clonalityTab[[3]], na.rm = TRUE)
	}
}

#output table
propTabName <- paste(subSample[1,6], holdingDir, "divergencyTab.txt", sep="")
#plotTab <- read.table(file=propTabName, sep="\t", header=TRUE)
#row.names(plotTab) <- c("trunk", "branch", "leaf")

#reorder list
clonalListCopy <- clonalList
reorderNum <- c((order(plotTab[1, c(12:16)]) + 11), (order(plotTab[1, c(20:23)]) + 19), (order(plotTab[1, c(17:19)]) + 16), order(plotTab[1,c(1:11)]), (order(plotTab[1, c(24:27)]) + 23)) 
counter <- 1
for(reOrd in reorderNum){
  clonalList[[counter]] <- clonalListCopy[[reOrd]]
  counter <- counter + 1
}

#order samples for plotting into adenoma / carcinoma groups
plotTab <- plotTab[reorderNum]
write.table(plotTab, file=propTabName, sep="\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
#reorderNum <- c(order(plotTab[1, c(1:9)]), 10:27)
#plotTab <- plotTab[reorderNum]

#plot clonal distributions
pdf(file=paste(subSample[1,6], plotDir, "divergencyPlot.pdf", sep=""), onefile=TRUE, width=6, height=4)
par(xpd=TRUE, mar=c(2,6,2,2))

#start plotting
plot(1, 1, col="white", axes=F, xlim=c(0,1), ylim=c(0,length(sampleNames)), xlab="", ylab="", main="")

#for each sample set plot rectangles
for(currPlot in 1:length(sampleNames)){
  #trunk rect
  rect(xleft = 0, xright = plotTab[1, currPlot], ybottom = (currPlot-1), ytop = currPlot, col="steelblue", border = "white")
  
  #branch rect
  rect(xleft = plotTab[1, currPlot], xright = (plotTab[1, currPlot] + plotTab[2, currPlot]), ybottom = (currPlot-1), ytop = currPlot, col="goldenrod", border = "white")
  
  #leaf rect
  rect(xleft = (plotTab[1, currPlot] + plotTab[2, currPlot]), xright = 1, ybottom = (currPlot-1), ytop = currPlot, col="salmon", border = "white")
  
  #add confidence lines if subsampling was performed
  tempDataHolder <- clonalList[[currPlot]]
  tempDataHolder <- tempDataHolder[complete.cases(tempDataHolder), ]
  if(nrow(tempDataHolder) > 1){
    upperLine <- quantile(tempDataHolder[[1]])["25%"]
    lowerLine <- quantile(tempDataHolder[[1]])["75%"]
    lines(x=c(upperLine, lowerLine), y=c((currPlot-0.5), (currPlot-0.5)), col="black")
    
    #upperBranLine <- min(tempDataHolder[[2]]) + plotTab[1, currPlot]
    #lowerBranLine <- max(tempDataHolder[[2]]) + plotTab[1, currPlot]
    #lines(x=c(upperBranLine, lowerBranLine), y=c((currPlot-0.35), (currPlot-0.35)))
    
    upperLeafLine <- quantile(tempDataHolder[[3]])["25%"] + plotTab[1, currPlot] + plotTab[2, currPlot]
    lowerLeafLine <- quantile(tempDataHolder[[3]])["75%"] + plotTab[1, currPlot] + plotTab[2, currPlot]
    lines(x=c(upperLeafLine, lowerLeafLine), y=c((currPlot-0.15), (currPlot-0.15)), col="black")
  }
  
  #add sample set names
  text(x=-0.01, y=(currPlot-0.5), labels = names(plotTab)[currPlot], pos = 2)
}

#add mean truncal
meanAdenoma <- mean(as.numeric(plotTab[1, c(1:5)]))
upperAdQrt <- quantile(as.numeric(plotTab[1, c(1:5)]), na.rm = TRUE)["75%"]
lowerAdQrt <- quantile(as.numeric(plotTab[1, c(1:5)]), na.rm = TRUE)["25%"]

#meanLynch <- mean(as.numeric(plotTab[1, c(24:27)]))
#lynchupQrt <- quantile(as.numeric(plotTab[1, c(24:27)]), na.rm = TRUE)["75%"]
#lynchlowQrt <- quantile(as.numeric(plotTab[1, c(24:27)]), na.rm = TRUE)["25%"]

meanCancer <- mean(as.numeric(plotTab[1, c(13:25,27)]))
upperQrt <- quantile(as.numeric(plotTab[1, c(13:25,27)]), na.rm = TRUE)["75%"]
lowerQrt <- quantile(as.numeric(plotTab[1, c(13:25,27)]), na.rm = TRUE)["25%"]


#plot mean lines
lines(x=c(meanAdenoma, meanAdenoma), y=c(1, 5), col="grey", lwd=2, lty=1)
lines(x=c(upperAdQrt, upperAdQrt), y=c(1, 5), col="red", lwd=2, lty=1)
lines(x=c(lowerAdQrt, lowerAdQrt), y=c(1, 5), col="red", lwd=2, lty=1)

#lines(x=c(meanLynch, meanLynch), y=c(24, 27), col="grey", lwd=2, lty=1)
#lines(x=c(lynchupQrt, lynchupQrt), y=c(24, 27), col="red", lwd=2, lty=1)
#lines(x=c(lynchlowQrt, lynchlowQrt), y=c(24, 27), col="red", lwd=2, lty=1)

lines(x=c(meanCancer, meanCancer), y=c(13, 23), col="grey", lwd=2, lty=1)
lines(x=c(upperQrt, upperQrt), y=c(13, 23), col="red", lwd=2, lty=1)
lines(x=c(lowerQrt, lowerQrt), y=c(13, 23), col="red", lwd=2, lty=1)

axis(side=1, at=seq(0,1,0.25), labels=seq(0,1,0.25))

dev.off()


#perform stats

#adenoma vs Lynch
#ks.test(as.numeric(plotTab[1, c(24:26)]), as.numeric(plotTab[1, c(1:5)]))

#adenoma vs carcinoma & Lynch (excluding 2 biopsy sets)
ks.test(as.numeric(plotTab[1, c(13:25,27)]), as.numeric(plotTab[1, c(1:5)]))

#carcinoma vs Lynch
#ks.test(as.numeric(plotTab[1, c(13:23)]), as.numeric(plotTab[1, c(24:26)]))
