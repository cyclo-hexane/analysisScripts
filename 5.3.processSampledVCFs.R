# 5.3.processSampledVCFs

################### notes: #####################

# processing order:
#
#
# 5.0.countBamReads [run on Apocrita, for each set] (gets tables of read numbers for a given bed file )
#        |
#        v
# 5.1.subsetBamFiles (calculates proportion to be subsampled for each set and makes .sh scripts to subset bams) 
#        |
#        v
# coverageBed run on samples for specified RIO as part of above script
#        |
#        v
# 5.2.sampleROIcoverage (assessed coverage then samples consistent regions 100 times for 10,000 regions)
#        |
#        v
# runPlatypus on sampled regions for subsetted bams
#        |
#        v
# 5.3.processSampledVCFs (to get somatic variants and count interval sampled by each iteration)
#        |
#        v
# 5.4.plotSampleDiversity
#
#################### libraries #################

annoVarList <- c()
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

#################### main program #################


#### main program ####

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
holdingDir <- "3.samplingAnalysis/platypusCalls/"
regionsDir <- "3.samplingAnalysis/sampledRegions/"
vcfName <- ".merged.vcf"

#concatenate names in table and delete unwanted columns
sampleList[2] <- paste(sampleList[[6]], holdingDir, sampleList[[1]],"/", sampleList[[1]], sep="")

#remove unwanted columns and make unique set table
sampleList <- sampleList[-c(3:4)]
sampleList <- unique(sampleList[1:6])

#filter each sampled vcf
for(j in 1:nrow(sampleList)){
  noSamples <- sampleList[j,6]
  
  for(currIter in 1:100){
    print(paste("#### filtering sample ", sampleList[j,1], " iteration ", currIter, " ####",sep=""))
    
    #### these file are deleted ####
    FILTName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".", currIter, ".FILT.vcf", sep="")
    VARName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".", currIter, ".VAR.vcf", sep="")
    SOMAName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".", currIter, ".SOMA.vcf", sep="")	
    logName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".", currIter, ".log", sep="")
    confName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".", currIter, ".somatic.vcf", sep="")
    
    #somatic files
    confOutput <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".", currIter, ".somatic.txt", sep="")
    confTotalName <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".", currIter, ".somatic.total.vcf", sep="")
    confTotalOutput <- paste(sampleList[j,4], holdingDir, sampleList[j,1],"/", sampleList[j,1], ".", currIter, ".somatic.total.txt", sep="")
    
    #prepare sample field strings
    counterTemp <- 1
    indexStrings <- as.list(NA)
    for(k in 0:(noSamples-1)){
      indexStrings[[counterTemp]] <- paste(" (GEN[", k,"].NR > 9) & ", sep="")
      counterTemp <- counterTemp +1
    }
    indexStrings[[length(indexStrings)]] <- substr(indexStrings[[length(indexStrings)]], 1,  (nchar(indexStrings[[length(indexStrings)]]) - 2 ))
    
    counterTemp <- 1
    snvStrings <- as.list(NA)
    for(k in 0:(noSamples-1)){
      snvStrings[[counterTemp]] <- paste(" (GEN[", k,"].NV > 0 ) | ", sep="")
      counterTemp <- counterTemp +1
    }
    snvStrings[[length(snvStrings)]] <- substr(snvStrings[[length(snvStrings)]], 1,  (nchar(snvStrings[[length(snvStrings)]]) - 2 ))
    
    #prepare indexes for extractFields command
    counterTemp <- 1
    extractStringsNR <- as.list(NA)
    for(k in 0:(noSamples-1)){
      extractStringsNR[[counterTemp]] <- paste(" \"GEN[", k,"].NR\" ", sep="")
      counterTemp <- counterTemp +1
      
    }
    
    counterTemp <- 1
    extractStringsNV <- as.list(NA)
    for(k in 0:(noSamples-1)){
      extractStringsNV[[counterTemp]] <- paste(" \"GEN[", k,"].NV\" ", sep="")
      counterTemp <- counterTemp +1
    }
    
    
    # 1 .filter by FILTER field
    filtVarCommand <- paste("cat ", sampleList[j,2], ".", currIter, vcfName, " | java -jar ~/bin/SnpSift.jar filter \"( ( (FILTER ='PASS') | (FILTER ='alleleBias') | (FILTER ='HapScore') | (FILTER ='SC') | (FILTER ='SC;alleleBias') | (FILTER ='HapScore;alleleBias') | (FILTER ='HapScore;SC') ) )\" > ", FILTName, sep="")
    system(command=filtVarCommand)
    
    # 2. annotate variant types
    annoCommand <- paste("java -jar bin/SnpSift.jar varType ", FILTName," > ", VARName, sep="")
    system(command=annoCommand)
    
    # 3. filter for somatic single nucleotide variants (not in normal)
    somaticCommand <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( ( (exists SNP) & (GEN[", sampleList[j,5],"].NV = 0) & (GEN[", sampleList[j,5],"].NR < 40) & (GEN[", sampleList[j,5],"].NR > 9) ) | ( (exists SNP) & (GEN[", sampleList[j,5],"].NV < 4) & (GEN[", sampleList[j,5],"].NR > 39) ) )\" > ", SOMAName, sep="")
    system(command=somaticCommand)
    
    # 4. filter by read depth (>9X for all samples), this is conformance somatic single nucleotide variant file
    depthCommand <- paste("cat ", SOMAName, " | java -jar ~/bin/SnpSift.jar filter \"( ( ", paste(indexStrings, collapse=" "), " ) & (", paste(snvStrings, collapse=" "), " ) )\" > ", confName, sep="")
    system(command=depthCommand)
  
    extractCommand <- paste("java -jar ~/bin/SnpSift.jar extractFields ", confName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", confOutput, sep="")
    system(command=extractCommand)
    
    #tidy up
    system(command=paste("rm ", FILTName, sep=""))
    system(command=paste("rm ", VARName, sep=""))
    system(command=paste("rm ", SOMAName, sep=""))
    system(command=paste("rm ", confName, sep=""))
    system(command=paste("rm ", logName, sep=""))
  
  }
}

#loop through each bed file and count intervals sampled (in Mb)
regionIntervals <- data.frame(matrix(NA, nrow=100, ncol=1))
names(regionIntervals) <- "size(Mb)"

for(currIt in 1:100){
  
  currRegFile <- paste(sampleList[1,4], regionsDir,"SeqCap.EZ.subsetted.", currIt, ".bed", sep="")
  itIn <- read.table(file=currRegFile, sep="\t", header = FALSE)
  
  #count total interval size
  totalIntervalTemp <- 0
  for(currRow in 1:nrow(itIn)){
    totalIntervalTemp <- totalIntervalTemp + (itIn[currRow, 3] - itIn[currRow, 2])
  }
  
  totalIntervalTemp <- totalIntervalTemp / 1000000
  
  #add to results table
  regionIntervals[currIt, 1] <- totalIntervalTemp
}

#output region sizes files
fileOut <- paste(sampleList[1,4], regionsDir,"SeqCap.EZ.totalSizeOfIntervals.txt", sep="")
write.table(regionIntervals, file=fileOut, sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
