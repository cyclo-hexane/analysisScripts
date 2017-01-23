# uses genetic and physical distance measures to produce linear model

######################## notes ########################
#
#
#
#
######################## libraries ########################


######################## subroutines ########################


######################## main program ########################


#input sampleList from commandline arguments
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
holdingDir <- "1.platypusCalls/exome/"
namePrepended <- ".snv.annoVar.exonic_variant_function.txt"
plotDir <- "10.finalFigures/"

sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[-c(7,12:17)]

#load driver gene list
driverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table.05-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)


######################## assess driver stats ########################


disDir <- "11.TumourImages/spatialLinearRegression/"
outDir <- "11.TumourImages/divData/"

listCounter <- 1
divList <- as.list(NA)

#calculate driver stats for cancers
for(j in 1:length(sampleNames)){
	print(paste("#### analyzing sample ", sampleNames[j], " ####",sep=""))
	
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
  
	for(currAF in 1:noSamples){
	  dataIn[ncol(dataIn)+1] <- dataIn[[8+noSamples+currAF]] / dataIn[[8+currAF]]
	}
	
	names(dataIn) <- c("line", "type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(samNames, ".NV", sep=""), paste(samNames, ".NR", sep=""), samNames)
  
  #get distance measures
  disIn <- read.csv(file=paste(subSample[1,6], disDir, setName, ".csv", sep=""), header=TRUE, stringsAsFactors = FALSE)
  names(disIn) <- c("samples", "Distance", "geneDiv", "totalVar")
  
  #for each pair get samples and assess divergence
  for(currRow in 1:nrow(disIn)){
    tempNames <- strsplit(disIn[currRow, 1], split=":")
    
    tempVars <- dataIn[c(tempNames[[1]])]
    tempVars[tempVars > 0] <- 1
    tempVars[ncol(tempVars)+1] <- tempVars[1] + tempVars[2]
    disIn[currRow, "totalVar"] <- nrow(tempVars[tempVars[[3]]!=0,])
    disIn[currRow, "geneDiv"] <- nrow(tempVars[tempVars[[3]]==1,])
  }
  
  disInTemp <- disIn
  disInTemp["maxDist"] <- max(disInTemp[["Distance"]])
  divList[[listCounter]] <- disInTemp
  listCounter <- listCounter + 1
  
  #output new tables
  #outFile <- paste(subSample[1,6], outDir, setName, ".divTab.csv", sep="")
  #write.csv(disIn, file=outFile, quote = FALSE, row.names = FALSE)
  
  #normalize values
  #disIn["Distance"] <- disIn[["Distance"]] / max(disInTemp[["Distance"]])
  #disIn["totalVar"] <- disIn[["totalVar"]] / max(disIn[["totalVar"]])
  
  #perform linear model
  #lmResults <- lm(data=disIn, formula = Distance ~ 0 + geneDiv)
  lmResults <- lm(data=disIn, formula = Distance ~ geneDiv)
  nd <- data.frame(x=seq(0, (2*max(disIn[["geneDiv"]])), length=100))
  names(nd) <- "geneDiv"
  pConf1 <- predict(lmResults, interval="prediction", newdata=nd)
  #pConf2 <- predict(lmResults, interval="confidence", newdata=nd)
  
  pdf(file=paste(subSample[1,6], "11.tumourImages/linearModel/", setName, ".divModel.pdf", sep=""), onefile=TRUE, width=5, height=5)
  par(mar=c(4,4,4,4))
  plot(x=disIn[["geneDiv"]], y=disIn[["Distance"]], pch=20, xlim=c(0,(2*max(disIn[["geneDiv"]]))), ylim=c(0,(2*max(disIn[["Distance"]]))), col="red", cex=2, main="physical distance vs genetic distance", ylab="physical dist", xlab="divergent exonic SNVs")
  
  abline(lmResults)
  text(x=(0.8*max(disIn[["geneDiv"]])), y=(1.8*max(disIn[["Distance"]])), labels = paste("R = ", round(summary(lmResults)$r.squared, digits = 5)), pos=4)
  text(x=(0.8*max(disIn[["geneDiv"]])), y=(1.7*max(disIn[["Distance"]])), labels = paste("p =", round(summary(lmResults)$coefficients[1,4], digits = 5)), pos=4)
  
  #confidence lines
  matlines(nd, pConf1[,c("lwr","upr")], col="grey50", lty=2, type="p", pch=1, cex=0.5)
  #matlines(nd, pConf2[,c("lwr","upr")], col="violet", lty=2, type="p", pch=1, cex=0.5)
  dev.off()
  
}





#perform linear regression
lmTab <- data.frame(matrix(NA, nrow=0, ncol=5))
names(lmTab) <- names(disIn)
for(currSam in c(1:3,5:length(divList))){
  lmTab <- rbind(lmTab, divList[[currSam]])
}
lmTab["Distance"] <- lmTab[["Distance"]] / lmTab[["maxDist"]] 
lmResults <- lm(data = lmTab, formula = Distance ~ 0 + geneDiv)
nd <- data.frame(x=seq(0, (2*max(disIn[["geneDiv"]])), length=100))
names(nd) <- "geneDiv"
pConf1 <- predict(lmResults, interval="prediction", newdata=nd)
#pConf2 <- predict(lmResults, interval="confidence", newdata=nd)

pdf(file=paste(subSample[1,6], "11.tumourImages/linearModel/total.divModel.pdf", sep=""), onefile=TRUE, width=5, height=5)
par(mar=c(4,4,4,4))
plot(x=lmTab[["geneDiv"]], y=lmTab[["Distance"]], pch=20, xlim=c(0,(2*max(lmTab[["geneDiv"]]))), ylim=c(0,(2*max(lmTab[["Distance"]]))), col="red", cex=1, main="physical distance vs genetic distance", ylab="physical dist", xlab="divergent exonic SNVs")

abline(lmResults)
text(x=(0.8*max(lmTab[["geneDiv"]])), y=(1.8*max(lmTab[["Distance"]])), labels = paste("R = ", round(summary(lmResults)$r.squared, digits = 5)), pos=4)
text(x=(0.8*max(lmTab[["geneDiv"]])), y=(1.7*max(lmTab[["Distance"]])), labels = paste("p =", round(summary(lmResults)$coefficients[1,4], digits = 5)), pos=4)

#confidence lines
matlines(nd, pConf1[,c("lwr","upr")], col="grey50", lty=2, type="p", pch=1, cex=0.2)
#matlines(nd, pConf2[,c("lwr","upr")], col="violet", lty=2, type="p", pch=1, cex=0.5)
dev.off()

