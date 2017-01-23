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
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.csv", header=FALSE, stringsAsFactors=FALSE)
holdingDir <- "1.platypusCalls/exome/"
namePrepended <- ".snv.annoVar.exonic_variant_function.txt"
plotDir <- "1.platypusCalls/exome/"

sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[c(18:20)]

#load driver gene list
driverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table5-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)


######################## assess SNV clonality ########################

#plot table
plotTab <- data.frame(matrix(NA, ncol=length(sampleNames), nrow=3))
names(plotTab) <- sampleNames
row.names(plotTab) <- c("trunk", "adenoma", "cancer")

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
	names(clonalityTab) <- c("trunk", "adenoma", "cancer")
  
  #count variants
  trunkCounter <- 0
  adenoma1 <- 0
  adenoma2 <- 0
  leafCounter <- 0
  for(currRow in 1:nrow(dataIn)){
    varTemp <- dataIn[currRow, samNamesNoNorm] > 0
    
    if(as.numeric(table(varTemp)["TRUE"]) == 3){
      trunkCounter <- trunkCounter + 1
    }else if(varTemp[1] == "FALSE" & varTemp[2] == "TRUE"){
      adenoma1 <- adenoma1 + 1
    }else if(varTemp[1] == "FALSE" & varTemp[3] == "TRUE"){
      adenoma2 <- adenoma2 + 1
    }else if(varTemp[1] == "TRUE" & varTemp[2] == "FALSE" & varTemp[3] == "FALSE"){
      leafCounter <- leafCounter + 1
    }
  } 
  #total variants
  branchCounter <- mean(adenoma1, adenoma2)
  totalVar <- trunkCounter + branchCounter + leafCounter
  
  #add proportions to clonality table
  clonalityTab[1, "trunk"] <- trunkCounter / totalVar
  clonalityTab[1, "adenoma"] <- branchCounter / totalVar
  clonalityTab[1, "cancer"] <- leafCounter / totalVar
  
  #store clonality table
  clonalList[[j]] <- clonalityTab
  
  #add averages to plot tab
  plotTab[1, j] <- mean(clonalityTab[[1]], na.rm = TRUE)
  plotTab[2, j] <- mean(clonalityTab[[2]], na.rm = TRUE)
  plotTab[3, j] <- mean(clonalityTab[[3]], na.rm = TRUE)
}


#order samples for plotting into adenoma / carcinoma groups
plotTab <- plotTab[reorderNum]

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

#add mean truncal and leaf lines
meanAdenoma <- mean(as.numeric(plotTab[1, c(12:16)]))
upperAdQrt <- quantile(as.numeric(plotTab[1, c(12:16)]), na.rm = TRUE)["75%"]
lowerAdQrt <- quantile(as.numeric(plotTab[1, c(12:16)]), na.rm = TRUE)["25%"]

meanCancer <- mean(as.numeric(plotTab[1, c(1:11)]))
upperQrt <- quantile(as.numeric(plotTab[1, c(1:11)]), na.rm = TRUE)["75%"]
lowerQrt <- quantile(as.numeric(plotTab[1, c(1:11)]), na.rm = TRUE)["25%"]

lines(x=c(meanCancer, meanCancer), y=c(1, 11), col="red", lwd=2, lty=2)
lines(x=c(upperQrt, upperQrt), y=c(1, 11), col="red", lwd=2, lty=2)
lines(x=c(lowerQrt, lowerQrt), y=c(1, 11), col="red", lwd=2, lty=2)

lines(x=c(meanAdenoma, meanAdenoma), y=c(11, 16), col="red", lwd=2, lty=2)
lines(x=c(upperAdQrt, upperAdQrt), y=c(11, 16), col="red", lwd=2, lty=2)
lines(x=c(lowerAdQrt, lowerAdQrt), y=c(11, 16), col="red", lwd=2, lty=2)

axis(side=1, at=seq(0,1,0.25), labels=seq(0,1,0.25))

dev.off()


