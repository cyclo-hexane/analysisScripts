# makes annovar input files and produces a heat map of variants
# takes input sample list as argument and assumes holding directory is output
# in its current state this script only plots SNVs, CNV infomration has been hashed out 
 
# ammended version

#### main program ####

require(gplots)
require(GMD)

#get latest version of heatmap
#source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#input sampleList from commandline arguments
arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=4){
	stop("\n#### please use syntax > Rscript 4B.variantHeatmaps.R < sample list file > < holding directory > < prepended.name > < homoplasy tables > ####\n")
}
sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
holdingDir <- arguments[2]
namePrepended <- arguments[3]
holdingDir2 <- arguments[4]

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
namePrepended <- ".snv.annoVar.exonic_variant_function.0.01"
holdingDir <- "1.platypusCalls/nonSyn.0.01/"
holdingDir2 <- "2.phylogenetics/nonSyn.0.01/"
sampleNames <- unique(sampleList[[1]])

#load driver gene list
driverList <- read.csv(file="~/PhD/CRCproject/10.finalFigures/supp.table5-driverRef.csv", header=TRUE, stringsAsFactors=FALSE)

#now process .vcf file and make table 
for(j in 1:length(sampleNames)){
	print(paste("#### making heatmap for sample ", sampleNames[j], " ####",sep=""))
	
	#subset main list
	subSample <- subset(sampleList, sampleList[1]==sampleNames[j])
	
	setName <- unique(subSample[[1]])
	noSamples <- subSample[1,8]
	normalIndex <- 1+subSample[1,7]
	
	
	#setup input/output names
	dataIn <- read.table(file=paste(subSample[1,6], holdingDir, setName,"/", setName, namePrepended,".txt", sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
	

	#filter out synonymous variants
	#dataIn <- subset(dataIn, dataIn[2]=="nonsynonymous SNV" | dataIn[2]=="stopgain")
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
	
	#make blank column for ordering value
	tempTable[ncol(tempTable)+1] <- NA
	tempTable[ncol(tempTable)+1] <- NA
	
	#order variants by clonal catagory placing value in end column and using sort
	for(k in 1:nrow(tempTable)){
		assessRow <- tempTable[k,9:(8+noSamples)] !=0
		assessTable <- table(assessRow)
		
		#catch bulk variants not in set 
		if(names(assessTable)[1]=="FALSE" & as.integer(assessTable["FALSE"])==noSamples){
			tempTable[k,(ncol(tempTable)-1)] <- NA
			next
		}
		#if trunkal mark as 3
		if(as.integer(assessTable["TRUE"])==noSamples){
			tempTable[k,(ncol(tempTable)-1)] <- 3
		}
		if(as.integer(assessTable["TRUE"])==1){
			tempTable[k,(ncol(tempTable)-1)] <- 1
		}
		if(as.integer(assessTable["TRUE"])!=1 & as.integer(assessTable["TRUE"])!=noSamples){
			tempTable[k,(ncol(tempTable)-1)] <- 2
		}
		
		#marks exact number of instances
		#tempTable[k,(ncol(tempTable)-1)] <- as.integer(assessTable["TRUE"])	
	}
	
	#remove bulk marked NA variants
	tempTable <- tempTable[!is.na(tempTable[ncol(tempTable)-1]), ]
	
	#now mark row AF sums for secondary ordering
	for(k in 1:nrow(tempTable)){
		assessRow <- tempTable[k,9:(8+noSamples)]				
		tempTable[k,ncol(tempTable)] <- sum(assessRow)	
	}
	
	#final sort (clonal catagory and AF sum only)
	sortTable <- tempTable[ order( -tempTable[(ncol(tempTable)-1)], -tempTable[ncol(tempTable)] ), ]
  
	#change catagories to colour names
	sortTable[sortTable[(ncol(sortTable)-1)]==1,(ncol(sortTable)-1)] <- "salmon"
	sortTable[ sortTable[(ncol(sortTable)-1)]==2,(ncol(sortTable)-1)] <- "goldenrod"
	sortTable[sortTable[(ncol(sortTable)-1)]==3,(ncol(sortTable)-1)] <- "steelblue"
	
	sortTable <- sortTable[-2]
	
	#now parse gene names into new column
	for(l in 1:nrow(sortTable)){
		geneName <- strsplit(sortTable[l,2], ":")
		sortTable[l,2] <- geneName[[1]][1]
		if(sortTable[l,2]=="stopgain"){
			sortTable[l,2] <- paste(sortTable[l,2], " (", sortTable[l,2], ") ", sep="")
		}
	}
	
	#### get homoplasy number state ####
	
	#get homoplasy file
	homopFile <- paste(subSample[1,6], holdingDir2, subSample[1,1], "/", subSample[1,1], ".homoAnno.txt", sep="")
	
	#merge position fields in main table
	sortTable[1] <- paste(sortTable[[3]], ":", sortTable[[4]], sep="")
	
	if(file.exists(homopFile)){
	  homoFileIn <- read.table(file=homopFile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
	  homoFileIn[5] <- paste(homoFileIn[[3]], ":", homoFileIn[[4]], sep="")
  	
  	#add to column 5 on sortTable, homoplasy matches, leave blank if not
  	for(m in 1:nrow(sortTable)){
  		if(sortTable[m, 1] %in% homoFileIn[[5]]){
  			sortTable[m, 5] <- "purple"
  		}else{
  			sortTable[m, 5] <- "white"
  		}	 
  	}
	}else{
	  sortTable[5] <- "white"
	}
		
		
	write.table(sortTable, file=paste(subSample[1,6], holdingDir, setName,"/", subSample[1,1],".exome.ordered2.txt", sep=""), sep="\t", quote=FALSE)
	
  
	#order for publication (by appearence in sample)
	
	#save trunk variants
	trunkVar <- sortTable[sortTable[(7+noSamples+1)]=="steelblue", ]
	branchVar <- sortTable[sortTable[(7+noSamples+1)]=="goldenrod" & sortTable[5]=="white", ]
	homoplasyVar <- sortTable[sortTable[(7+noSamples+1)]=="goldenrod" & sortTable[5]=="purple", ]
	leafOrder <- sortTable[sortTable[(7+noSamples+1)]=="salmon", ]
  
	samCol <- c(8:(7+noSamples))
  
  #order branch and homoplasy variants
	orderTemp <- branchVar[c(8:(7+noSamples))]
	if(nrow(orderTemp)!=0){
  	orderTemp[orderTemp > 0] <- 1
  	orderTemp[ncol(orderTemp)+1] <- NA
  	for(currOrder in 1:nrow(orderTemp)){
  	  orderTemp[currOrder, ncol(orderTemp)] <- paste(orderTemp[currOrder, c(1:noSamples)], collapse = ".")
  	}
  	branchVar <- branchVar[order(orderTemp[ncol(orderTemp)], decreasing = TRUE),]
	}
	
	orderTemp <- homoplasyVar[c(8:(7+noSamples))]
	if(nrow(orderTemp)!=0){
  	orderTemp[orderTemp > 0] <- 1
  	orderTemp[ncol(orderTemp)+1] <- NA
  	for(currOrder in 1:nrow(orderTemp)){
  	  orderTemp[currOrder, ncol(orderTemp)] <- paste(orderTemp[currOrder, c(1:noSamples)], collapse = ".")
  	}
  	homoplasyVar <- homoplasyVar[order(orderTemp[ncol(orderTemp)], decreasing = TRUE),]
  	homoplasyVar[(7+noSamples+1)] <- "grey"
	}
  
	#order leaf variants
  for(currOrder in samCol){
    leafOrder <- leafOrder[order(leafOrder[currOrder], decreasing = TRUE), ]
  }
  
  #re-assemble variant table
	sortTable <- rbind(trunkVar, branchVar, homoplasyVar, leafOrder)
  
	#get driver gene color scheme for plot
	driverCols <- data.frame(matrix(NA, ncol=1, nrow=nrow(sortTable)))
	driverCols[1] <- "black"
	for(n in 1:nrow(sortTable)){
	  geneTemp <- strsplit(sortTable[n,2], ":")
	  if(geneTemp[[1]][1] %in% driverList[[1]] | geneTemp[[1]][1] %in% driverList[[2]]){
	    driverCols[n,1] <- "red"
	  }
	}
	
	
	#### plot via rectangle ####
	
	outputLoc <- paste(subSample[1,6], holdingDir, subSample[1,1], "/", subSample[1,1], ".AFmap.v2.pdf", sep="")
	pdf(file=outputLoc, width=(0.09*nrow(tempTable)), height=noSamples+1)
	par(xpd=TRUE, mar=c(4,5,8,3))
	plot(1, 1, col="white", axes=F, xlim=c(0,nrow(sortTable)), ylim=c(0,noSamples*1.9),xlab="", ylab="", main="")
	
	#get sample name list
	rectNameList <- subSample[-normalIndex,2]
	#rectNameList <- rev(rectNameList)
	
	ycount <- 1.9
	normalCounter <- 1
	for(cancer in 8:(7+(noSamples))){
				
		#reset gene counter
		xcount <- 1
		
		#get current sample AF data and colours
		currentSample <- sortTable[c(2,cancer, (ncol(sortTable)-1))]
		currentSample[ncol(currentSample)+1] <- sortTable[5]
		
		#get copy number table for this sample
		#cnvIn <- read.table(file=paste(subSample[1,6], "CNVs/snv.gt.calls/", rectNameList[normalCounter],".cloneHD.force1.wl0.1.penalty0.95.snv.gt.txt.snv.gt.txt", sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE, colClasses=rep("numeric",4) )
		
		#change chromosome names for X and Y and concatenate name
		# cnvIn[cnvIn[[1]]==23, 1] <- "X"
		# cnvIn[cnvIn[[1]]==24, 1] <- "Y"
		# cnvIn[1] <- paste(cnvIn[[1]], ":", cnvIn[[2]], sep="")
		
		# #loop thought and assign states
		# cnvColList <- data.frame(matrix(0, nrow=nrow(sortTable), ncol=3))
		# for(getCNV in 1:nrow(sortTable)){
			# if(currentSample[getCNV,4] %in% cnvIn[[1]]){
				# cnvColList[getCNV,1] <- (cnvIn[cnvIn[[1]]==currentSample[getCNV,4] ,3])*0.2
				# cnvColList[getCNV,2] <- (cnvIn[cnvIn[[1]]==currentSample[getCNV,4] ,5])*0.2
				
				# #set maximum states for plot
				# if(cnvColList[getCNV,2] > 0.8){
					# cnvColList[getCNV,2] <- 0.8
					# if(cnvColList[getCNV,1] > 0.8){
						# cnvColList[getCNV,1] <- 0.8
					# }
				# }
				 
				# if(cnvColList[getCNV,2] <= 0.3 & cnvColList[getCNV,2] > 0){
					# cnvColList[getCNV,3] <- "red"
				# }
				# if(cnvColList[getCNV,2] >= 0.3 & cnvColList[getCNV,2] < 0.55 & cnvColList[getCNV,1] <= 0.3 & cnvColList[getCNV,1] <= 0.3){
					# cnvColList[getCNV,3] <- "grey75"
				# }
				# if(cnvColList[getCNV,2] >= 0.55 & cnvColList[getCNV,2] < 0.75){
					# cnvColList[getCNV,3] <- "blue"
				# }
				# if(cnvColList[getCNV,2] >= 0.75){
					# cnvColList[getCNV,3] <- "blue4"
				# }
				# if(cnvColList[getCNV,2] <= 0.45 & cnvColList[getCNV,1] <= 0.45 & cnvColList[getCNV,2] >= 0.3 & cnvColList[getCNV,1] >= 0.3){
					# cnvColList[getCNV,3] <- "red"
				# }		
			# }else{
				# cnvColList[getCNV,3] <- "white"
			# } 
		# }
		
		for(gene in 1:nrow(currentSample)){
			
			#get allele frequency of current variant
			currentAF <- currentSample[gene,2] 
			if(currentAF>0 & currentAF<0.85){
				rect(xleft=xcount ,xright=xcount+1 ,ybottom=ycount ,ytop=(ycount+currentAF) ,col=currentSample[gene,3], border ="white")
				rect(xleft=xcount ,xright=xcount+1 ,ybottom=(ycount+currentAF) ,ytop=(ycount+1) ,col=paste("light", currentSample[gene,3], sep=""), border ="white")
			}
			if(currentAF>0.84){
				rect(xleft=xcount ,xright=xcount+1 ,ybottom=ycount ,ytop=(ycount+currentAF) ,col=paste(currentSample[gene,3],4,sep=""), border = "white")
				rect(xleft=xcount ,xright=xcount+1 ,ybottom=(ycount+currentAF) ,ytop=(ycount+1) ,col=paste("light", currentSample[gene,3], sep=""), border = "white")
			}
			if(currentAF==0){
				rect(xleft=xcount ,xright=xcount+1 ,ybottom=ycount ,ytop=(ycount+currentAF) ,col="white", border = "white")
			}
				
			#plot copy number state section	
			#rect(xleft=xcount ,xright=xcount+1 ,ybottom=ycount-0.8 ,ytop=(ycount-0.79+cnvColList[gene,2]) ,col=cnvColList[gene,3], border ="white")
			#rect(xleft=xcount ,xright=xcount+1 ,ybottom=ycount-0.8 ,ytop=(ycount-0.79+cnvColList[gene,1]) ,col="grey50", border ="white")
			
			xcount <- xcount + 1
		}
		#mark dotter and paritioning lines between graphs
		lines(x=c(1,nrow(currentSample)+1) ,y=c(ycount+0.5, ycount+0.5), lty=2)
		lines(x=c(1,nrow(currentSample)+1) ,y=c(ycount+1, ycount+1), lty=1)
		lines(x=c(1,nrow(currentSample)+1) ,y=c(ycount-0.01, ycount-0.01), lty=1)
		lines(x=c(1,nrow(currentSample)+1) ,y=c(ycount-0.8-0.01, ycount-0.8-0.01), lty=1, lwd=0.5)
		lines(x=c(1,nrow(currentSample)+1) ,y=c(ycount-0.8-0.01, ycount-0.8-0.01), lty=1, lwd=0.5)
		lines(x=c(1,nrow(currentSample)+1) ,y=c(ycount-0.8+0.4, ycount-0.8+0.4), lty=2, lwd=0.5)
		lines(x=c(1,nrow(currentSample)+1) ,y=c(ycount-0.8+0.6, ycount-0.8+0.6), lty=2, lwd=0.5)
		lines(x=c(1,nrow(currentSample)+1) ,y=c(ycount-0.8+0.2, ycount-0.8+0.2), lty=2, lwd=0.5)
		
		#label up sample names
		axis(2 ,at= ycount+0.5, labels=rectNameList[normalCounter] ,lwd=0 ,las=2 ,cex.axis=1 ,line=-1.5)
		
		ycount <- ycount + 1.9
		normalCounter <- normalCounter+1
	}
	
	#make homoplasy marks
	for(homoplasyMark in 1:nrow(sortTable)){
		rect(xleft=homoplasyMark, xright=homoplasyMark+1, ybottom=(ycount-0.85), ytop=(ycount-0.7), col=sortTable[homoplasyMark,5], border ="white")
	}
	
	#label up gene names in loop (to get different colors)
	for(p in 1:nrow(currentSample)){
		axis(3 ,at=p+0.5 ,labels=currentSample[[1]][p] ,lwd=0 ,las=2 ,cex.axis=0.4, line=0.5, col.axis=driverCols[[1]][p])
	}
	
	#label up allele frequency axis
	axis(4 ,at=ycount-1.9+0.5, labels=0.5 ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-1.8)
	axis(4 ,at=ycount-1.9+0.94, labels=1 ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-1.8)
	axis(4 ,at=ycount-1.9+0.1, labels=0 ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-1.8)
	axis(4 ,at=ycount-1.9+0.5, labels="AF" ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-0.95)
	
	#label up copy number bargraph axis
	axis(4 ,at=ycount-1.9-0.8+0.1, labels=1 ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-1.8)
	axis(4 ,at=ycount-1.9-0.8+0.3, labels=2 ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-1.8)
	axis(4 ,at=ycount-1.9-0.8+0.5, labels=3 ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-1.8)
	axis(4 ,at=ycount-1.9-0.8+0.7, labels=4 ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-1.8)
	axis(4 ,at=ycount-1.9-0.8+0.4, labels="CN state" ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-0.95)

	
	#label homoplasy row
	axis(4 ,at=ycount-0.73, labels="homoplasy" ,lwd=0 ,las=2 ,cex.axis=0.5 ,line=-0.95)
	
	#reset names for publication
	#rectNameList <- c("Left Crypt 1" ,"Left Crypt 2" ,"Left Crypt 3" ,"Left Crypt 4" ,"Left Crypt 5" ,"Right Crypt 1" ,"Right Crypt 2" ,"Right Crypt 3" ,"Right Crypt 4" ,"Right Crypt 5" ,"Right Crypt 6")
	
	dev.off()
}