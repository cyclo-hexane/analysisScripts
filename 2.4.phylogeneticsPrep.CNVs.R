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




################## main program ################## 


sampleList <- read.csv(file="~/PhD/CRCproject/8.CNVphylogenetics/trainingSampleList.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/archive/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)


holdingDir <- "8.CNVphylogenetics/"

vcfName <- ".phyloCNVs.csv"

apocritaDir <- "phylogenetics/CNVphylo/"

setNames <- unique(sampleList[[1]])

#loop through samples and get binary strings
for(x in 1:length(setNames)){
	
	#get sample names
	currentNames <- subset(sampleList, sampleList[1]==setNames[x])
	sampleNames <- currentNames[[2]]
  
	print(paste("making tree building scripts for " , currentNames[1,1] ,sep=""))
  
	noSamples <- currentNames[1,8]
	normalIndexTab <- currentNames[1,7]+1
	sampleNames <- sampleNames[-normalIndexTab]
	normalIndex <- noSamples
  
	CNVfileIn <- paste(currentNames[1,6], holdingDir, currentNames[1,1], vcfName, sep="")
	confData <- read.table(file=CNVfileIn, sep=",", header=TRUE, stringsAsFactors=FALSE)
	
  #remove those that cannot be modelled
	confData <- confData[confData[["model"]]!="cannotModel", ]
  
  copyConf <- confData
	
	#remove unwanted columns and convert AF to binary	
  confData <- confData[sampleNames]
	confData <- as.data.frame(append(confData, 0, after = currentNames[1,7]))
  names(confData)[normalIndexTab] <- currentNames[normalIndexTab, 2]
	
	#save genotypes as binary strings in list
	sampleData <- as.list(0)
	for(z in 1:ncol(confData)){		
		sampleData[[z]] <- confData[[z]]
	}
	
	#get output locations apocrita .tre log files and then .nex file local output respectively
	outputLoc <- paste(currentNames[1,5], apocritaDir, currentNames[1,1], "/", currentNames[1,1], sep="")
	outputLocNex <- paste(currentNames[1,6], holdingDir, currentNames[1,1], "/", currentNames[1,1], sep="")
	histogramFile <- paste(currentNames[1,5], apocritaDir, currentNames[1,1], "/", currentNames[1,1], ".hist.txt", sep="")
	
	#store strings with names in list
	counter<-0
	asStrings <- as.list(0)
	for(k in seq(1,(2*noSamples),2)){
		counter <- counter+1
		asStrings[[k]] <- as.character(currentNames[counter, 2])
		asStrings[[k+1]] <- paste(sampleData[[counter]], collapse='')
	}
	
	#### create file names
	noTaxa <- noSamples
	taxaNames <- paste(currentNames[[2]], collapse=" ")
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

	outgroup ", currentNames[normalIndexTab,2]," /only;

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

	outgroup ", currentNames[normalIndexTab,2]," /only;

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

	outgroup ", currentNames[normalIndexTab,2]," /only;

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