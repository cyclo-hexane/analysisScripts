# script processes .vcf (converted using snpsift to .txt) ready for phylogenetics
# 1. gets sample set file locations from argument and indivdual sample names from name changes file
# 2. convert allele frequencies into binary sequence for nexus file
# 3. outputs two nexus files (tree generation and bootsrapping), and apocrita .sh script

# updated version

################## notes #######################
# .txt (.vcf) Files listed in above have be in format:
# CHROM, POS, REF, ALT, GEN[x].NR, GEN[x].NV, 4 x Annotation fields.
# changes list is the same file used in program 1.processVCFsnpSift.R
# only used to get sample names (line 51) Example: ~/PhD/CRCproject/sampleLists/nameChangesList.total.txt


################## subroutines ################## 

filterVAF <- function(confDataSub, filtSub, normalIndexSub, outFileNameSub){
	#get rows to remove
	removeData <- confDataSub	
	removeData[normalIndexSub] <- 0
	samples <- c(1:noSamples)
	samples <- samples[-normalIndexSub]
	
	#assign phylogenetic states
	for(currRow in 1:nrow(removeData)){
		tempTab <- table(removeData[currRow,]>0)
		if(as.integer(tempTab["TRUE"]) != (noSamples-1) & as.integer(tempTab["TRUE"]) != 1){
			for(currSam in samples){
				if(removeData[currRow, currSam] <= filtSub){
					removeData[currRow, currSam] <- 0
				}
			}
		}
	}
	
	#check for zero columns and mark to remove
	for(currRow in 1:nrow(removeData)){
		if(sum(removeData[currRow,]) == 0){
			removeData[currRow, normalIndexSub] <- 1
		}
	}
	
	#filter sets appropriately
	removeVector <- which(removeData[normalIndexSub] > 0)
	if(length(removeVector) > 0 ){
		keepVector <- c(1:nrow(removeData))
		keepVector <- keepVector[-removeVector]
		keepVector <- removeData[keepVector,]
		
		#output new file
		write.table(removeData, file=outFileNameSub, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
		return(keepVector)
	}else{
		#output new file
		write.table(removeData, file=outFileNameSub, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
		return(removeData)
	}
}

filterVAFbyCNV <- function(confDataSub, filtSub, normalIndexSub, noSamplesSub){
	#convert X and Y chromosomes to 23,24
	chrCol <- as.character(confDataSub[[4]])
	chrCol[chrCol == "X"] <- 23
	confDataSub[4] <- chrCol
	chrCol <- as.character(confDataSub[[4]])
	chrCol[chrCol == "Y"] <- 24
	confDataSub[4] <- chrCol
	confDataSub <- subset(confDataSub, confDataSub[4]==1 | confDataSub[4]==2 | confDataSub[4]==3 | confDataSub[4]==4 | confDataSub[4]==5 | confDataSub[4]==6 | confDataSub[4]==7 | confDataSub[4]==8 | confDataSub[4]==9 | confDataSub[4]==10 | confDataSub[4]==11 | confDataSub[4]==12 | confDataSub[4]==13 | confDataSub[4]==14 | confDataSub[4]==15 | confDataSub[4]==16 | confDataSub[4]==17 | confDataSub[4]==18 | confDataSub[4]==19 | confDataSub[4]==20 | confDataSub[4]==21 | confDataSub[4]==22 | confDataSub[4]==23 | confDataSub[4]==24 )
	
	#get cnv file
	cnvData <- read.table(file=filtSub, sep="\t", header=TRUE, stringsAsFactors=FALSE)	
	
	#mark unknowns as wild type
	cnvData[cnvData == "X"] <- 1
	
	#remove wild type intervals
	cnvRemove <- c(0)
	counterSub <- 1
	for(currRow in 1:nrow(cnvData)){
		tempTable <- table(cnvData[currRow, c(5:(4 + (noSamplesSub-1) * 2))] == 1)
		if("TRUE" %in% names(tempTable)){
			if(tempTable["TRUE"] == ((noSamplesSub-1) * 2)){
				cnvRemove[counterSub] <- currRow
				counterSub <- counterSub + 1
			}
		}
		
	}
	if(cnvRemove[1] != 0){
		cnvData <- cnvData[-cnvRemove,]
	}
	
	#get SNV row to remove
	snvRemove <- c()
	counterSub <- 1
	for(currRow in 1:nrow(confDataSub)){
		if(as.numeric(confDataSub[currRow,4]) %in% cnvData[[1]]){
			subCNV <- subset(cnvData, cnvData[1]==confDataSub[currRow,4])	
		}else{
			next
		}
		for(currCNV in 1:nrow(subCNV)){
			if(confDataSub[currRow,5] > subCNV[currCNV,2] & confDataSub[currRow,5] < subCNV[currCNV,4]){
				snvRemove[counterSub] <- currRow
				counterSub <- counterSub + 1
			}
		}
	}
	snvRemove <- unique(snvRemove)
  #confDataSub <- confDataSub[-snvRemove,]
	return(snvRemove)	
}



################## main program ################## 



arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=5){
	stop("\\usage: > Rscript 2.phylogeneticsPrep.annoVar.R <fileList.csv> <holding directory> <output dir> <prependedName.vcf> <apocrita holding directory>")
}

sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.Set.10.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- arguments[2]
holdingDir2 <- arguments[3]
#holdingDir <- "1.platypusCalls/somaticTotal/"
#holdingDir2 <- "2.phylogenetics/somaticTotal/"

vcfName <- arguments[4]
#vcfName <- ".snv.annoVar.exonic_variant_function.txt"
#vcfName <- ".snv.annoVar.variant_function.txt"

#vcfOut <- ".filteredAF.txt"

#add filter vector
#filterVector <- c(0.08, 0.04, 0.04, 0.04, 0.06, 0.10, 0.80, 0.00, 0.04, 0.08, 0.06, 0.10, 0.10, 0.10, 0.10, 0.04)

apocritaDir <- arguments[5]
#apocritaDir <- "phylogenetics/somaticTotal.0.01/"
#apocritaDir <- ""

#get samplenames list
sampleNames <- sampleList[1:2]

#concatenate names in table and delete unwanted columns
sampleList[2] <- paste(sampleList[[6]], holdingDir, sampleList[[1]],"/", sampleList[[1]], vcfName, sep="")

#remove unwanted columns and make unique set table
sampleList <- sampleList[-c(3:4)]
sampleList <- unique(sampleList[1:6])

typeFlag <- "genome"

#loop through samples and get binary strings
for(x in 1:nrow(sampleList)){
	#input .txt file
	print(paste("making tree building scripts for sample " ,sampleList[x,1] ,sep=""))
	
	#get sample names
	currentNames <- subset(sampleNames, sampleNames[1]==sampleList[x,1], select=V2)
	
	noSamples <- sampleList[x,6]
	normalIndex <- sampleList[x,5]+1
	
	confData <- read.table(file=sampleList[x,2], sep="\t", header=FALSE, stringsAsFactors=FALSE)
	
	if(typeFlag == "genome"){
	  confData <- as.data.frame(append(confData, NA, after=0))
	}
  
  copyConf <- confData
  
	#calculate allele frequencies
	for(y in 1:noSamples){
		confData[(ncol(confData)+1)] <- as.vector(confData[(8+noSamples+y)]) / as.vector(confData[8+y])
	}
	
	#remove unwanted columns and convert AF to binary	
	confData <- confData[,-c(9:(8+noSamples*2))]
	
	#filter for non-synonymous mutations (if exome)
	#confData <- subset(confData, confData[2]=="nonsynonymous SNV" | confData[2]=="stopgain")
	
	#filter for exonic mutations (if genome)
	#confData <- subset(confData, confData[2]=="exonic")
	
	#filter by CNV regions
	#CNVfile <- paste("~/PhD/CRCprojectBackup/CNVs/cloneHD.2015.05.01/cloneHD-copynumber-states-v1.0/", sampleList[x,1], ".penalty0.90.baf.gt.txt", sep="")
	#filtVar <- filterVAFbyCNV(confData, CNVfile, normalIndex, noSamples)
	
	#write table for plotting reference
	#refOut <- paste(sampleList[x,4], holdingDir2, sampleList[x,1],"/", sampleList[x,1],".noCNV.filt.0.05.txt", sep="")
	#noCNVtab <- copyConf[-filtVar,]
	#confData <- confData[-filtVar,]
  #write.table(noCNVtab, file=refOut, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)	
	
	#remove unwanted column
	confData <- confData[-c(1:8)]
	
	#filter data by allele frequency
	#outFileName <- paste(sampleList[x,4], holdingDir, sampleList[x,1],"/", sampleList[x,1], vcfOut, sep="")
	#confData <- filterVAF(confData, filterVector[x], normalIndex, outFileName)
	
	#zero can be changed here for any inclusion criteria for phylogenetics, such as AF > 0.1 etc
	confData[confData > 0] <- 1
	#confData[confData > filterVector[x]] <- 1
	#confData[confData <= filterVector[x]] <- 0
	
	indexLoc <- as.integer(sampleList[x,5])
	indexLoc <- indexLoc+1
	
	#zero off normal column (hash out if AFs to be retained)
	confData[indexLoc] <- 0
	
	#save genotypes as binary strings in list
	sampleData <- as.list(0)
	for(z in 1:ncol(confData)){		
		sampleData[[z]] <- confData[[z]]
	}
	
	#get output locations apocrita .tre log files and then .nex file local output respectively
	outputLoc <- paste(sampleList[x,3], apocritaDir, sampleList[x,1], "/", sampleList[x,1], sep="")
	outputLocNex <- paste(sampleList[x,4], holdingDir2, sampleList[x,1], "/", sampleList[x,1], sep="")
	histogramFile <- paste(sampleList[x,3], apocritaDir, sampleList[x,1], "/", sampleList[x,1], ".hist.txt", sep="")
	
	#store strings with names in list
	counter<-0
	asStrings <- as.list(0)
	for(k in seq(1,(2*noSamples),2)){
		counter <- counter+1
		asStrings[[k]] <- as.character(currentNames[counter,1])
		asStrings[[k+1]] <- paste(sampleData[[counter]], collapse='')
	}
	
	#### create file names
	noTaxa <- noSamples
	taxaNames <- paste(currentNames[[1]], collapse=" ")
	noCharacters <- nchar(asStrings[[2]][[1]])
	matrixStrings <- paste(asStrings, collapse="\n")
	ascestState <- 0
	upperHIvalue <- round(noCharacters + ((noCharacters / 100) * 55), 0)
	for(w in 1:(noCharacters-1)){ascestState <- paste(ascestState,0,sep="")}
	treeFile <- paste(outputLoc, ".tre", sep="")
	treeFileAll <- paste(outputLoc, ".allTreeSearch.tre", sep="")
	logFile <- paste(outputLoc, ".log", sep="")
	logFileAll <- paste(outputLoc, ".allTrees.log", sep="")
	logFileBoot <- paste(outputLoc, ".boot.log", sep="")
	nexusTreeFileLoc <- paste(outputLoc, ".nex", sep="")
	nexusBootFileLoc <- paste(outputLoc, ".boot.nex", sep="")
	nexusAllTreesFile <- paste(outputLoc, ".allTrees.nex", sep="")

	
	#put all strings into nexus format for tree production
	totalStringsTree <- paste("#NEXUS
begin paup;
  set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
  log start file= ", logFile,"  replace;

	begin taxa;
        dimensions ntax=", noTaxa,";
        taxlabels
       ", taxaNames,";
	end;
	
	begin characters;
        dimensions nchar= ", noCharacters," ;
        format symbols = \"01\";
        matrix\n",matrixStrings,";

	end;
	
	begin assumptions;
        options deftype = unord;

        usertype statetransitionsmatrix (stepmatrix)= 2
		0 1 
		. i 
		1 . 
		;

	end;
		
	hsearch addseq=simple nreps=20 hold=1 rearrlimit=10000000 nbest=1000;

	outgroup ", currentNames[indexLoc,1]," /only;

	roottrees outroot=monophyl rootmethod=outgroup;

	describetrees /plot=both root=outgroup brlens=yes chglist=yes apolist=yes diag=yes homoplasy=yes xout=both;

	savetrees from=1 to=1000 file= ", treeFile,"  savebootp=Both format=altnex brlens=yes root=yes append=yes;
  
log stop;

end;
		
quit;", sep="")
	
	
	#write nexus file
	lapply(totalStringsTree, write, paste(outputLocNex,".nex",sep=""), append=FALSE)
	
	#put all strings into nexus format for all tree production
	totalStringsTreeAll <- paste("#NEXUS
begin paup;
  set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
  log start file= ", logFileAll,"  replace;

	begin taxa;
        dimensions ntax=", noTaxa,";
        taxlabels
       ", taxaNames,";
	end;
	
	begin characters;
        dimensions nchar= ", noCharacters," ;
        format symbols = \"01\";
        matrix\n",matrixStrings,";

	end;
	
	begin assumptions;
        options deftype = unord;

        usertype statetransitionsmatrix (stepmatrix)= 2
		0 1 
		. i 
		1 . 
		;

	end;
		
	alltrees fdfile=", histogramFile," keep=", upperHIvalue,";

	outgroup ", currentNames[indexLoc,1]," /only;

	roottrees outroot=monophyl rootmethod=outgroup;

	describetrees /plot=both root=outgroup brlens=yes chglist=yes apolist=yes diag=yes homoplasy=yes xout=both;

	savetrees from=1 to=1000 file= ", treeFileAll,"  savebootp=Both format=altnex brlens=yes root=yes append=yes;
  
log stop;

end;
		
quit;", sep="")
	
	
	#write nexus file
	lapply(totalStringsTreeAll, write, paste(outputLocNex,".allTrees.nex",sep=""), append=FALSE)

	
	
	#put all strings into nexus format for bootstrap production
	totalStringsBoot <- paste("#NEXUS
begin paup;
  set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
  log start file= ", logFileBoot,"  replace;

	begin taxa;
        dimensions ntax=", noTaxa,";
        taxlabels
       ", taxaNames,";
	end;
	
	begin characters;
        dimensions nchar= ", noCharacters," ;
        format symbols = \"01\";
        matrix\n",matrixStrings,";

	end;
	
	begin assumptions;
        options deftype = unord;

        usertype statetransitionsmatrix (stepmatrix)= 2
		0 1 
		. i 
		1 . 
		;

	end;
		
	hsearch addseq=simple nreps=20 hold=1 rearrlimit=10000000;

	outgroup ", currentNames[indexLoc,1]," /only;

	roottrees outroot=monophyl rootmethod=outgroup;

	bootstrap nreps=10000;

	log stop;

end;
		
quit;", sep="")
	
	
	#write nexus file
	lapply(totalStringsBoot, write, paste(outputLocNex,".boot.nex",sep=""), append=FALSE)

	#create shell scripts to run nexus files
	shellStrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 2            # Request 2 CPU cores
#$ -l h_rt=48:0:0      # Request 48 hour runtime
#$ -l h_vmem=11G    # Request 12GB RAM / core, i.e. 24GB total

./bin/paup4b10-x86-linux -n ", nexusTreeFileLoc, "
./bin/paup4b10-x86-linux -n ",nexusBootFileLoc, "
./bin/paup4b10-x86-linux -n ", nexusAllTreesFile
, sep="")

	#write nexus file
	lapply(shellStrings, write, paste(outputLocNex,".runPAUP.sh",sep=""), append=FALSE)

}
