# 8.0.msiSensorScripts

################### notes: #####################


#################### libraries #################


#################### main program #################

#sample list
sampleList <- read.csv(file="~/PhD/CRCproject/3.samplingAnalysis/masterSampleList.sampling.csv", header=FALSE, stringsAsFactors=FALSE)

#prepended names
prepName <- "_statistics.sample_interval_summary"
holdingDir <- "3.samplingAnalysis/coverageFiles/"
apocritaDir <- "processedBams/subSample/sampledRegions"

#remove unwanted samples
sampleList <- sampleList[sampleList[[1]]!="Polyp.08", ]

setNames <- unique(sampleList[[1]])

#make msi sensor scripts
for(currSam in 1:length(setNames)){ 
  #get normal file name
  

  MSIOutFile <- paste(sampleList[currSet, 6], shellDir, "runPlatypus.", currSamName, ".", currSam,".sh", sep="")
    
  totalStringsSH <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 6            # Request 6 CPU cores
#$ -l h_rt=48:0:0      # Request 48 hour runtime
#$ -l h_vmem=3.5G    # Request 4GB RAM / core, i.e. 24GB total
                          
msisensor msi \
-d /data/BCI-EvoCa/william/referenceHG19/hs37d5.microsalellites.list \
-n /data/ \
-t /data/ \
-o /data/ \
-b 6" 
, sep="")
    
  lapply(totalStringsSH, write, MSIOutFile, append=FALSE)
  
}


