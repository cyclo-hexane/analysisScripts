# script processes .vcf (converted using snpsift to .txt) ready for phylogenetics
# in this version the allele frequency cutoff is stratified for each sample
# 1. gets sample set file locations from argument and indivdual sample names from name changes file
# 2. convert allele frequencies into binary sequence for nexus file
# 3. outputs two nexus files (tree generation and bootsrapping), and apocrita .sh script

# updated version

#### notes ####
# .txt (.vcf) Files listed in above have be in format:
# CHROM, POS, REF, ALT, GEN[x].NR, GEN[x].NV, 4 x Annotation fields.
# changes list is the same file used in program 1.processVCFsnpSift.R
# only used to get sample names (line 51) Example: ~/PhD/CRCproject/sampleLists/nameChangesList.total.txt


arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=4){
	stop("\\usage: > Rscript 2.phylogeneticsPrep.annoVar.strat.R < fileList.csv > < holding directory > < .prependedName.vcf > < apocrita holding directory >")
}

sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- arguments[2]
#holdingDir <- "1.platypusCalls/rawData/"

vcfName <- arguments[3]
#vcfName <- ".snv.annoVar.exonic_variant_function.txt"

apocritaDir <- arguments[4]
#apocritaDir <- "2.phylogenetics/stratified/"

#get samplenames list
sampleNames <- sampleList[1:2]

#concatenate names in table and delete unwanted columns
sampleList[2] <- paste(sampleList[[6]], holdingDir, sampleList[[1]],"/", sampleList[[1]],vcfName, sep="")

#remove unwanted columns and make unique set table
sampleList <- sampleList[-c(3:4)]
sampleList <- unique(sampleList[1:6])

#loop through samples and get binary strings
for(x in 1:nrow(sampleList)){
	#input .txt file
	print(paste("making scripts for sample " ,sampleList[x,1] ,sep=""))
	
	newDir <- paste(sampleList[x,4], apocritaDir, sampleList[x,1], "/", sep="")
	system(command=paste("mkdir", newDir))
	
	dataIn <- read.table(file=sampleList[x,2], sep="\t", header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
	
	#get sample names
	currentNames <- subset(sampleNames, sampleNames[1]==sampleList[x,1], select=V2)

	noSamples <- sampleList[x,6]
	
	nexFileList <- c(0)
	nonBootList <- c(0)
  
	#filter for non-synonymous mutations
	dataIn <- subset(dataIn, dataIn[2]=="nonsynonymous SNV" | dataIn[2]=="stopgain")
	
	#switch for exome or full annotations (extra column for exome)
	dataIn <- dataIn[-c(1:3,5)]
	#dataIn <- dataIn[-c(1:2,5)]
	
	#calculate allele frequencies
	for(y in 1:noSamples){
		dataIn[(ncol(dataIn)+1)] <- dataIn[(4+noSamples+y)]/dataIn[4+y]
	}
	
	#remove unwanted columns and convert AF to binary	
	dataIn <- dataIn[,-c(1:(4+(noSamples*2)))]
	
  #for each AF filter, generate 100 bootstrap trees
	counter2 <- 1
  counter3 <- 1
	for(w in seq(0,0.35,0.02)){
		print(paste("making script with filter " ,w ,sep=""))
		
		confData <- dataIn
		
		#filter out by allele frequency
		confData[confData > w] <- 1
		confData[confData <= w] <- 0
		
		#remove columns containing all zeros
		confData <- confData[apply(confData, 1, function(x) !all(x==0)),]
		
		indexLoc <- as.integer(sampleList[x,5])
		indexLoc <- indexLoc+1
		
		#zero off normal column (hash out if AFs to be retained)
		confData[indexLoc] <- 0
		
		sampleData <- as.list(0)
		for(z in 1:ncol(confData)){  	
		  sampleData[[z]] <- confData[[z]]
		}
    
		#get output locations apocrita .tre log files and then .nex file local output respectively
		outputLoc <- paste(sampleList[x,3], apocritaDir, sampleList[x,1], "/", sampleList[x,1], "_f", w, sep="")
		outputLocNex <- paste(sampleList[x,4], apocritaDir, sampleList[x,1], "/", sampleList[x,1], "_f", w, sep="")
		
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
		treeFile <- paste(outputLoc, ".tre", sep="")
		logFile <- paste(outputLoc, ".log", sep="")
		nexusTreeFileLoc <- paste(outputLoc, ".nex", sep="")
		
		#put all strings into nexus format for tree production
		totalStringsTree <- paste("#NEXUS
begin paup;
set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no;
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
			
		hsearch addseq=simple nreps=20 hold=1 rearrlimit=10000000 nbest=10;
	
		outgroup ", currentNames[indexLoc,1]," /only;
	
		roottrees outroot=monophyl rootmethod=outgroup;
	
		describetrees /plot=both root=outgroup brlens=yes chglist=yes apolist=yes diag=yes homoplasy=yes xout=both;
	
		savetrees from=1 to=10 file= ", treeFile,"  savebootp=Both format=altnex brlens=yes root=yes append=yes;
	  
	log stop;
	
	end;
			
	quit;", sep="")
		
		
		#write nexus file
		lapply(totalStringsTree, write, paste(outputLocNex,".nex",sep=""), append=FALSE)
		
    nonBootList[counter3] <- paste(outputLoc, sep="")
		counter3 <- counter3+1
    
    
    ###### get each bootstrap file ######
    for(currBoot in 1:100){
      
      #get random selection of variants
      bootChar <- sample(c(1:nrow(confData)), nrow(confData), replace = TRUE)
      bootConf <- confData[bootChar, ]
      
      #save genotypes as binary strings in list
      sampleData <- as.list(0)
      for(z in 1:ncol(confData)){		
        sampleData[[z]] <- bootConf[[z]]
      }
      
      #get output locations apocrita .tre log files and then .nex file local output respectively
      outputLoc <- paste(sampleList[x,3], apocritaDir, sampleList[x,1], "/", sampleList[x,1], "_f", w, "_", currBoot, sep="")
      outputLocNex <- paste(sampleList[x,4], apocritaDir, sampleList[x,1], "/", sampleList[x,1], "_f", w, "_", currBoot, sep="")
      
      #save file name
      nexFileList[counter2] <- outputLoc
      counter2 <- counter2+1
      
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
      treeFile <- paste(outputLoc, ".tre", sep="")
      logFile <- paste(outputLoc, ".log", sep="")
      nexusTreeFileLoc <- paste(outputLoc, ".nex", sep="")
      
      #put all strings into nexus format for tree production
      totalStringsTree <- paste("#NEXUS
	begin paup;
	  set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no;
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
			
		hsearch addseq=simple nreps=20 hold=1 rearrlimit=10000000 nbest=10;
	
		outgroup ", currentNames[indexLoc,1]," /only;
	
		roottrees outroot=monophyl rootmethod=outgroup;
	
		describetrees /plot=both root=outgroup brlens=yes chglist=yes apolist=yes diag=yes homoplasy=yes xout=both;
	
		savetrees from=1 to=10 file= ", treeFile,"  savebootp=Both format=altnex brlens=yes root=yes append=yes;
	  
	log stop;
	
	end;
			
	quit;", sep="")
      
      
      #write nexus file
      lapply(totalStringsTree, write, paste(outputLocNex,".nex",sep=""), append=FALSE) 
    } 
	}
		
	#create shell scripts to run nexus files
	nexListRunFiles <- paste("./bin/paup4b10-x86-linux ", nexFileList, ".nex", sep="", collapse="\n")
	nexListBootFiles <- paste("./bin/paup4b10-x86-linux ", nexFileList, ".boot.nex", sep="", collapse="\n")
	outputLocNex <- paste(sampleList[x,4], apocritaDir, "/", sampleList[x,1], sep="")
	
	shellStrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 2            # Request 2 CPU cores
#$ -l h_rt=48:0:0      # Request 48 hour runtime
#$ -l h_vmem=12G    # Request 12GB RAM / core, i.e. 24GB total
	
", nexListRunFiles , sep="")
	
	#write nexus file
	lapply(shellStrings, write, paste(outputLocNex,".runPAUP.sh",sep=""), append=FALSE)
  
  
  
	#create shell scripts to run non-bootstrap nexus files
	nexListRunFiles <- paste("./bin/paup4b10-x86-linux ", nonBootList, ".nex", sep="", collapse="\n")
	outputLocNex <- paste(sampleList[x,4], apocritaDir, sampleList[x,1], sep="")
	
	shellStrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 2            # Request 2 CPU cores
#$ -l h_rt=48:0:0      # Request 48 hour runtime
#$ -l h_vmem=12G    # Request 12GB RAM / core, i.e. 24GB total
	
", nexListRunFiles , sep="")
	
	#write nexus file
	lapply(shellStrings, write, paste(outputLocNex,".nonBoot.runPAUP.sh",sep=""), append=FALSE)	
  
}

