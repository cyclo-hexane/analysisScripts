# script processes a list of .vcf files
# subroutine to rename sample genotype fields in .vcf not in use
# performs the following
# 1. filters by FILTER field i.e PASS/FAIL etc
# 2. annotate variant type (SNP/DEL/INS)
# 3. filter for somatic variants, not in normal (SNP and INDEL handled separately)
# 4. by read depth (>9 X for all samples) and not for total variants (non-conformance output)
# 5. annotate variant effects using snpEFF (hash out to ignore)
# 6. output .txt version of annotate muation .vcf

# updated version

#### notes ####
# uses standardized file list. Example: ~/PhD/cryptProject/masterCryptList.csv
# in subroutine gsed command is the GNU version on MacOSX which works properly! If using on Linux change to sed
# assumes multi-allele variants have been dealt with prior to running

#### subroutines ####

makeAnnovar <- function(samples, prepend, output){
	sampleNames <- unique(samples[[1]])
	namePrepended <- prepend

	#now process .vcf file and make annoVar input 
	for(j in 1:length(sampleNames)){
		print(paste("#### making table for sample ", sampleNames[j], " ####",sep=""))
	
		#subset main list
		subSample <- subset(samples, samples[1]==sampleNames[j])
		#print(subSample)
	
		setName <- unique(subSample[[1]])
	
		#setup input/output names
		dataIn <- read.table(file=paste(subSample[1,6], holdingDir, setName,"/", setName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE)
	
		dataIn <- as.data.frame(append(dataIn, list(to = dataIn[[2]]), after = 2))
	
		annoInput <- paste(subSample[1,6], holdingDir, setName,"/", setName, output,".annoVarInput.txt", sep="")
	
		write.table(dataIn, file= annoInput, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	
		annoOut <- paste(subSample[1,6], holdingDir, setName,"/", setName, output,".annoVar", sep="") 
	
		#run annovar program on new file
		system(command=paste("~/bin/annovar/annotate_variation.pl -out " ,annoOut ," -build hg19 ", annoInput, " ~/bin/annovar/humandb/",sep=""))
	
		#annovar output files
		annoExoFunct <- paste(subSample[1,6], holdingDir, setName,"/", setName, output,".annoVar.exonic_variant_function", sep="")
		annoExoFunctNew <- paste(subSample[1,6], holdingDir, setName,"/", setName, output,".annoVar.exonic_variant_function.txt", sep="")

		annoVarFunct <- paste(subSample[1,6], holdingDir, setName,"/", setName, output,".annoVar.variant_function", sep="")
		annoVarFunctNew <- paste(subSample[1,6], holdingDir, setName,"/", setName, output,".annoVar.variant_function.txt", sep="")

	
		#rename files with .txt prepend
		system(command=paste("mv ", annoExoFunct, " ", annoExoFunctNew, sep=""))
		system(command=paste("mv ", annoVarFunct, " ", annoVarFunctNew, sep=""))
	}
}

#### main program ####

arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
if(length(arguments)!=4){
	stop("\n#### please use syntax > Rscript 3.subsetVariants.R <sample list file> <holding directory> <prepended.name> <exclusions list> ####\n")
}

sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterLynchList.WGS.csv", header=FALSE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/PhD/CRCproject/masterAdinCaList.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- arguments[2]
#holdingDir <- "1.platypusCalls/rawData/"

vcfName <- arguments[3]
#vcfName <- ".merged.vcf"
#vcfName <- ".vcf"

sampleListCopy <- sampleList

#concatenate names in table and delete unwanted columns
sampleList[2] <- paste(sampleList[[6]], holdingDir, sampleList[[1]],"/", sampleList[[1]], vcfName, sep="")

#remove unwanted columns and make unique set table
sampleList <- sampleList[-c(3:4)]
sampleList <- unique(sampleList[1:6])

#get sample exclusion file
excList <- arguments[4]
#excList <- read.csv(file="~/PhD/CRCproject/masterLynchList.WGS.exc.csv", header=FALSE, stringsAsFactors=FALSE)
excList <- read.csv(file="~/PhD/CRCproject/masterAdinCaList.exc.csv", header=FALSE, stringsAsFactors=FALSE)

#perform pipeline for each sample
for(j in 1:nrow(sampleList)){
	print(paste("#### filtering sample ", sampleList[j,1], " ####",sep=""))
	
	#setup input/output names
	
	#### these file are deleted ####
	FILTName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".FILT.vcf", sep="")
	VARName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".VAR.vcf", sep="")
	SOMAName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".SOMA.vcf", sep="")
	indelSOMAName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".indel.SOMA.vcf", sep="")
	
	#### these file are kept for analysis ####
	villefilterName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".villeFilt.vcf", sep="")
	villefilterOutput <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".villeFilt.txt", sep="")
  
	#somatic files
	confName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".somatic.vcf", sep="")
	confOutput <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".somatic.txt", sep="")
	confTotalName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".somatic.total.vcf", sep="")
	confTotalOutput <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".somatic.total.txt", sep="")
	indelName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".indel.somatic.vcf", sep="")
	indelOutput <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".indel.somatic.txt", sep="")
	
	#SNP files
	snpName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".SNPs.vcf", sep="")
	snpOutput <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".SNPs.txt", sep="")
	snpTotalName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".SNPs.total.vcf", sep="")
	snpTotalOutput <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".SNPs.total.txt", sep="")
	
	#subset exclusion list
	excluSub <- subset(excList, excList[1]==sampleList[j,1])
	keepList <- (which(excluSub[7]==1))-1
	
	#remove normal column from keep list
	keepListNoNorm <- keepList[-which(keepList==sampleList[j,5])]
	
	#prepare indexes for depth and somatic filtering commands
	counterTemp <- 1
	indexStrings <- as.list(NA)
	for(k in keepListNoNorm){
		indexStrings[[counterTemp]] <- paste(" (GEN[", k,"].NR > 9) & ", sep="")
		counterTemp <- counterTemp +1
	}
	indexStrings[[length(indexStrings)]] <- substr(indexStrings[[length(indexStrings)]], 1,  (nchar(indexStrings[[length(indexStrings)]]) - 2 ))

	
	#prepare indexes for total variant output
	counterTemp <- 1
	totalIndexStrings <- as.list(NA)
	for(k in keepListNoNorm){
		totalIndexStrings[[counterTemp]] <- paste(" ((GEN[", k,"].NR > 9) & (GEN[", k,"].NV > 0)) | ", sep="")
		counterTemp <- counterTemp +1
	}
	totalIndexStrings[[length(totalIndexStrings)]] <- substr(totalIndexStrings[[length(totalIndexStrings)]], 1,  (nchar(totalIndexStrings[[length(totalIndexStrings)]]) - 2 ))

	
	counterTemp <- 1
	snvStrings <- as.list(NA)
	for(k in keepListNoNorm){
		snvStrings[[counterTemp]] <- paste(" (GEN[", k,"].NV > 0 ) | ", sep="")
		counterTemp <- counterTemp +1
	}
	snvStrings[[length(snvStrings)]] <- substr(snvStrings[[length(snvStrings)]], 1,  (nchar(snvStrings[[length(snvStrings)]]) - 2 ))
	
	counterTemp <- 1
	totalIndexSNPs <- as.list(NA)
	for(k in keepList){
		totalIndexSNPs[[counterTemp]] <- paste(" (GEN[", k,"].NR > 9) & ", sep="")
		counterTemp <- counterTemp +1

	}
	totalIndexSNPs[[length(totalIndexSNPs)]] <- substr(totalIndexSNPs[[length(totalIndexSNPs)]], 1,  (nchar(totalIndexSNPs[[length(totalIndexSNPs)]]) - 2 ))

	
	#prepare indexes for extractFields command
	counterTemp <- 1
	extractStringsNR <- as.list(NA)
	for(k in keepList){
		extractStringsNR[[counterTemp]] <- paste(" \"GEN[", k,"].NR\" ", sep="")
		counterTemp <- counterTemp +1

	}
	
	counterTemp <- 1
	extractStringsNV <- as.list(NA)
	for(k in keepList){
		extractStringsNV[[counterTemp]] <- paste(" \"GEN[", k,"].NV\" ", sep="")
		counterTemp <- counterTemp +1
	}
	

	# 1 .filter by FILTER field
	filtVarCommand <- paste("cat ", sampleList[j,2], " | java -jar ~/bin/SnpSift.jar filter \"( ( (FILTER ='PASS') | (FILTER ='alleleBias') | (FILTER ='HapScore') | (FILTER ='SC') | (FILTER ='SC;alleleBias') | (FILTER ='HapScore;alleleBias') | (FILTER ='HapScore;SC') ) )\" > ", FILTName, sep="")
	system(command=filtVarCommand)
	
	# 2. annotate variant types
	annoCommand <- paste("java -jar bin/SnpSift.jar varType ", FILTName," > ", VARName, sep="")
	system(command=annoCommand)
  
  
  #baf SNP files for ville
	#filterCommand <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( (exists SNP) & (FILTER ='PASS') & (isHet(GEN[", sampleList[j,5],"])) & (GEN[", sampleList[j,5], "].NR>20)  )\" > ", villefilterName, sep="")
	#system(command=filterCommand)
	
	#extractCommand6 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", villefilterName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", villefilterOutput, sep="")
	#system(command=extractCommand6)
	
	
	
	#### somatic files ####
	
	# 3. filter for somatic single nucleotide variants (not in normal)
	somaticCommand <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( ( (exists SNP) & (GEN[", sampleList[j,5],"].NV = 0) & (GEN[", sampleList[j,5],"].NR < 100) & (GEN[", sampleList[j,5],"].NR > 9) ) | ( (exists SNP) & (GEN[", sampleList[j,5],"].NV < 3) & (GEN[", sampleList[j,5],"].NR > 99) ) )\" > ", SOMAName, sep="")
	system(command=somaticCommand)
	
	# 3b. filter for somatic indels (not in normal)
	#somaticCommand2 <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( ( (exists INS) & (GEN[", sampleList[j,5],"].NV = 0) ) | ( (exists DEL) & (GEN[", sampleList[j,5],"].NV = 0) ) )\" > ", indelSOMAName, sep="")
	somaticCommand2 <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( exists INS ) | exists DEL )\" > ", indelSOMAName, sep="")
	system(command=somaticCommand2)

	# 4. filter by read depth (>9X for all samples), this is conformance somatic single nucleotide variant file
	depthCommand <- paste("cat ", SOMAName, " | java -jar ~/bin/SnpSift.jar filter \"( ( ", paste(indexStrings, collapse=" "), " ) & (", paste(snvStrings, collapse=" "), " ) )\" > ", confName, sep="")
	system(command=depthCommand)
	
	# 4. filter by read depth (>9X for ANY samples), this is total variants file
	depthCommand5 <- paste("cat ", SOMAName, " | java -jar ~/bin/SnpSift.jar filter \"( ( ", paste(totalIndexStrings, collapse=" ") ,"))\" > ", confTotalName, sep="")
	system(command=depthCommand5)
	
	# 4b. filter by read depth (>9X for all samples), this is conformance file INDELS
	depthCommand2 <- paste("cat ", indelSOMAName, " | java -jar ~/bin/SnpSift.jar filter \"( ( ", paste(indexStrings, collapse=" "), " ) & (", paste(snvStrings, collapse=" "), " ) )\" > ", indelName, sep="")
	system(command=depthCommand2)

	#### SNP files ####

	# 4c. produce non-fitered SNP file for total variant calculations
	depthCommand3 <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( (exists SNP) & (GEN[", sampleList[j,5],"].NV > 0) & (GEN[", sampleList[j,5],"].NR > 9) )\" > ", snpTotalName, sep="")
	system(command=depthCommand3)
	
	# 4d. produce fitered SNP file, this is SNP conformance
	depthCommand4 <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( (exists SNP) & ( ", paste(totalIndexSNPs, collapse=" ") ,") & (GEN[", sampleList[j,5],"].NV > 0) )\" > ", snpName, sep="")
	system(command=depthCommand4)


	#### vcf to text file conversion ####
	
	extractCommand <- paste("java -jar ~/bin/SnpSift.jar extractFields ", confName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", confOutput, sep="")
	system(command=extractCommand)
	
	extractCommand2 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", confTotalName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", confTotalOutput, sep="")
	system(command=extractCommand2)
		
	extractCommand3 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", indelName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", indelOutput, sep="")
	system(command=extractCommand3)
	
	extractCommand4 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", snpName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", snpOutput, sep="")
	system(command=extractCommand4)
	
	extractCommand5 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", snpTotalName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", snpTotalOutput, sep="")
	system(command=extractCommand5)

	#tidy up
	system(command=paste("rm ", FILTName, sep=""))
	system(command=paste("rm ", VARName, sep=""))
	system(command=paste("rm ", SOMAName, sep=""))
	system(command=paste("rm ", indelSOMAName, sep=""))

}

#annotate somatic file with annoVar
makeAnnovar(sampleListCopy, ".somatic.txt", ".snv")
makeAnnovar(sampleListCopy, ".somatic.total.txt", ".snv.total")
makeAnnovar(sampleListCopy, ".indel.somatic.txt", ".indel")
makeAnnovar(sampleListCopy, ".SNPs.txt", ".snps")
makeAnnovar(sampleListCopy, ".SNPs.total.txt", ".snps.total")


