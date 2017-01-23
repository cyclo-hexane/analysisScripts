# script processes annoVar file and filters by allele frequency


################## notes #######################
# 1. gets sample set file locations from argument and indivdual sample names from name changes file
# 2. filter exome file for nonsynon or synon variants (set at line 61-62)
# 3. filters by allele frequency (only branch variants, specified by line 48)
# 4. outputs files to holdingDir2

################ subroutines ###################


################ main program ###################

arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=4){
	stop("\\usage: > Rscript 1.2.processVCF.AFfilter < fileList.csv > < holding directory > < .prependedName.vcf > < vcf output name >")
}

sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/PhD/CRCproject/masterLynchList.WGS.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.Ha.01.filt.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- arguments[2]
#holdingDir <- "1.platypusCalls/somaticTotal/"
#holdingDir2 <- "1.platypusCalls/somaticTotal.0.01/"

vcfName <- arguments[3]
typeFlag <- "exome"
#typeFlag <- "genome"
#vcfName <- ".snv.annoVar.exonic_variant_function.txt"
#vcfName <- ".snv.annoVar.variant_function.txt"

vcfOut <- arguments[4]
#vcfOut <- ".snv.annoVar.exonic_variant_function.0.01.txt"
#vcfOut <- ".snv.annoVar.variant_function.0.01.txt"

#get samplenames list
setNames <- unique(sampleList[[1]])

#filter allele frequencies
filterVector <- rep(0.01, nrow(sampleList))
#filterVector <- 0.01

#depth cutoff for C->T variants
depthCO <- 20
#filterCutCT <- filterVector[1]
filterCutCT <- 0.1

#loop through samples and filter
for(x in 1:length(setNames)){
	print(paste("filtering sample " ,setNames[x] ,sep=""))
	
  subSample <- sampleList[sampleList[[1]]==setNames[x], ]
  sampleNames <- subSample[[2]]
  
	confData <- read.table(file=paste(subSample[1,6], holdingDir, setNames[x], "/", setNames[x], vcfName, sep=""), sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
	noSamples <- subSample[x, 8]
	normalIndex <- (sampleList[x, 7])+1
	
  if(typeFlag == "genome"){
    confData <- as.data.frame(append(confData, NA, after=0))
  }else{
    #blank off for full genome file
    #confData <- subset(confData, confData[2]=="nonsynonymous SNV" | confData[2]=="stopgain" | confData[2]=="stoploss")
    #confData <- subset(confData, confData[2]=="synonymous SNV")
  }

	#remove non standard chromosomes
	confData[4] <- as.character(confData[[4]])
	confData[7] <- as.character(confData[[7]])
	confData[8] <- as.character(confData[[8]])
	confData[confData[[4]]=="X", 4] <- "23"
	confData[confData[[4]]=="Y", 4] <- "24"
	#confData <- confData[confData[[4]] %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"), ]

	#calculate allele frequencies
	for(y in 1:noSamples){
		confData[(ncol(confData)+1)] <- confData[(8+noSamples+y)]/confData[8+y]
	}
	names(confData) <- c("line", "region", "variant", "chrom", "pos", "pos2", "ref", "alt", paste(subSample[[2]], ".NR", sep=""), paste(subSample[[2]], ".NV", sep=""), paste(subSample[[2]], ".AF", sep=""))
	
	#get rows to remove
	removeData <- c()
	remCounter <- 1
	samples <- subSample[[2]]
	samples <- samples[-normalIndex]
	
	####### filter by VAF filter cutoff, (and minimum read depth) 
	for(currRow in 1:nrow(confData)){
		tempRow <- confData[currRow, ]
		NVvalues <- tempRow[1, paste(samples, ".NV", sep="")]
		NRvalues <- tempRow[1, paste(samples, ".NR", sep="")]
		VAFvalues <- tempRow[1, paste(samples, ".AF", sep="")]
		currBases <- tempRow[1, c("ref", "alt")]
		
		#test read depth against curoff if a C->T variant
		if((as.character(currBases[1])=="C" & as.character(currBases[2])=="T") | (as.character(currBases[1])=="G" & as.character(currBases[2])=="A")){
  		
		  #test depths
		  depthTab <- table(NRvalues < depthCO)
		  
		  #if at least one biopsy has a low depth, test VAFs and depth individually
  		if("TRUE" %in% names(depthTab)){
  		  samCheckList <- rep(NA, length(samples))
  		  for(currTest in 1:length(samples)){
  		    if(as.numeric(VAFvalues[currTest]) > 0){
  		      #if sample has variant check depth and VAF
  		      if(as.numeric(VAFvalues[currTest]) > filterCutCT & as.numeric(NRvalues[currTest]) > depthCO){
  		        samCheckList[currTest] <- 1
  		      }else{
  		        samCheckList[currTest] <- 0
  		      }
  		    }else{
  		      #depth >10X and no variant present
  		      samCheckList[currTest] <- 1
  		    }
  		  }
  		  if(0 %in% samCheckList){
  		    #remove row
  		    removeData[remCounter] <- currRow
  		    remCounter <- remCounter + 1
  		    next()
  		  }
  		}else{
  		  removeData[remCounter] <- currRow
  		  remCounter <- remCounter + 1
  		  next()
  		}
		}else{
		  #test if clonal
		  clonalAssess <- table(VAFvalues > 0)
		  if(as.numeric(clonalAssess["TRUE"]) != (noSamples-1)){
		    #subclonal, assess each sample in turn and change to zero if below cutoff
		    for(currTest in 1:length(samples)){
		      if(as.numeric(VAFvalues[currTest]) < filterVector[x]){
		        #VAF below cutoff
		        confData[currRow, names(VAFvalues[currTest])] <- 0
		      }
		    }
		    
		    #mark to remove if all zeros
		    if(sum(confData[currRow, names(VAFvalues)]) == 0){
		      #remove row
		      removeData[remCounter] <- currRow
		      remCounter <- remCounter + 1
		      next()
		    }
		  }else{
		    #clonal, assess whether only one sample is below VAF cutoff
		    VAFassess <- table(VAFvalues > filterVector[x])
		    if("TRUE" %in% names(VAFassess)){
		      #assess whether low VAF in onyl one sample
		      if(as.numeric(VAFassess["TRUE"]) == (noSamples-1)){
		        #true clonal, do nothing
		      }else{
		        for(currTest in 1:length(samples)){
		          if(as.numeric(VAFvalues[currTest]) < filterVector[x]){
		            #VAF below cutoff
		            confData[currRow, names(VAFvalues[currTest])] <- 0
		          }
		        }
		        
		        #mark to remove if all zeros
		        if(sum(confData[currRow, names(VAFvalues)]) == 0){
		          #remove row
		          removeData[remCounter] <- currRow
		          remCounter <- remCounter + 1
		          next()
		        }
		      }
		    }else{
		      #subclonal in all samples, check indiviudally
		      for(currTest in 1:length(samples)){
		        if(as.numeric(VAFvalues[currTest]) < filterVector[x]){
		          #VAF below cutoff
		          confData[currRow, names(VAFvalues[currTest])] <- 0
		        }
		      }
		      
		      #mark to remove if all zeros
		      if(sum(confData[currRow, names(VAFvalues)]) == 0){
		        #remove row
		        removeData[remCounter] <- currRow
		        remCounter <- remCounter + 1
		        next()
		      }
		    } 
		  }
		}
	}
	
	
	#filter sets appropriately
	if(length(removeData) > 0){
		confRem <- confData[removeData, ]
		keepConfData <- confData[-removeData, ]
		
		#output new file
		outFileName <- paste(subSample[1,6], holdingDir2, subSample[1,1],"/", subSample[1,1], vcfOut, sep="")
		write.table(keepConfData["line", "region", "variant", "chrom", "pos", "pos2", "ref", "alt", paste(subSample[[2]], ".NR", sep=""), paste(subSample[[2]], ".NV", sep="")], file=outFileName, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
	
		#output removed data
		outFileRem <- paste(subSample[1,6], holdingDir2, subSample[1,1],"/", subSample[1,1], ".removed.txt", sep="")
		write.table(confRem["line", "region", "variant", "chrom", "pos", "pos2", "ref", "alt", paste(subSample[[2]], ".NR", sep=""), paste(subSample[[2]], ".NV", sep="")], file=outFileRem, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
	}else{
		#output new file
		outFileName <- paste(subSample[1,6], holdingDir2, subSample[1,1],"/", subSample[1,1], vcfOut, sep="")
		write.table(copyData["line", "region", "variant", "chrom", "pos", "pos2", "ref", "alt", paste(subSample[[2]], ".NR", sep=""), paste(subSample[[2]], ".NV", sep="")], file=outFileName, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
	}
}
