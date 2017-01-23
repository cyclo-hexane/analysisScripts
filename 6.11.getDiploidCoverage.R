#
#
##########################   notes   ##########################
#
#   1. reads segmentations tables and get diploid regions
#             |
#             V
#   2. makes bed file of these regions in 1Mb windows
#             |
#             V
#   3. outputs script to run on apocrita
#
#
########################## libraries  ##########################

######################### subroutines ##########################

######################### main program ##########################

#get sample list
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[1]!="Polyp.08", ]

setNames <- unique(sampleList[[1]])
setNames <- setNames[-c(1,4,6:8,12:14,16)]

#working dir
CNVholdingDir <- "7.CNVcalls.archive/baf.updated/"
CNVfileNames <- ".penalty0.95.baf.gt.txt"

#out directory
dirOut <- "7.CNVcalls.archive/runScripts/"


############ for each sample set ###############

for(currSam in 1:length(setNames)){
  #subset main list
  subSample <- subset(sampleList, sampleList[1]==setNames[currSam])
  
  sampleNames <- subSample[[2]]
  sampleNames <- sampleNames[-(subSample[1,7]+1)]
  noSamples <- length(sampleNames)
  
  #get seg file
  dataIn <- read.table(file=paste(subSample[1,6], CNVholdingDir, subSample[1,1], CNVfileNames, sep=""), sep="\t", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
  
  #get seg interval size
  dataIn[3] <- (dataIn[["last.locus"]] - dataIn[["first.locus"]]) / 1000000

  
  #assess each seg and assign Normal, Gain, Loss, counts and phylogenetic location
  for(currAss in 1:noSamples){
    print(paste("###### processing sample set ", subSample[1,1], " ", sampleNames[currAss], " ########", sep=""))
    
    majCol <- paste(sampleNames[currAss], "_Major", sep="")
    minCol <- paste(sampleNames[currAss], "_Minor", sep="")
    
    #make bed file for this sample
    bedName <- paste(subSample[1,6], dirOut, subSample[1,1], "_", subSample[currAss,2], ".diploidRegions.bed", sep="")
    bedFile <- data.frame(matrix(NA, nrow = 1000000, ncol=3))
    bedCounter <- 1
    
    for(assessRow in 1:nrow(dataIn)){
      currRow <- dataIn[assessRow, c(majCol, minCol)]
      if(currRow[1, majCol] == 1 & currRow[1, minCol] == 1){
        #add diploid segs to bed file
        assRow <- dataIn[assessRow, c(1:4)]
        if(assRow[1, "nloci"] >= 1){
          diploidBins <- seq(assRow[1,2], assRow[1,4], 1000000)
          if(length(diploidBins) > 1){
            for(addRow in 2:length(diploidBins)){
              bedFile[bedCounter, 1] <- assRow[1, "chr"] 
              bedFile[bedCounter, 2] <- diploidBins[(addRow-1)]
              bedFile[bedCounter, 3] <- diploidBins[addRow]-1
              bedCounter <- bedCounter + 1
            }
          }else{
            bedFile[bedCounter, 1] <- assRow[1, "chr"] 
            bedFile[bedCounter, 2] <- diploidBins[1]
            bedFile[bedCounter, 3] <- diploidBins[2]
            bedCounter <- bedCounter + 1
          }
        }else{
          #region too small to bin, exclude
          next
        }
      }else{
        #non-diploid seg, do nothing
      }
    }
    
    #output bed file
    bedFile <- bedFile[complete.cases(bedFile), ]
    write.table(bedFile, file=bedName, sep="\t", row.names = FALSE, col.names = FALSE)
    bamFile <- paste(subSample[1,5], "processedBams/", subSample[1,1], "/", subSample[currAss,2], ".mkdub.bam", sep="")
    outFile <- paste(subSample[1,5], "processedBams/", subSample[1,1], "/", subSample[currAss,2], "_diploid_statistics", sep="")
    bedNameApo <- paste(subSample[1,5], "processedBams/", subSample[1,1], "/", subSample[1,1], "_", subSample[currAss,2], ".diploidRegions.bed", sep="")
    
    #make target coverage file
    scriptName <- paste(subSample[1,6], dirOut, "runGATKcoverageBed.", subSample[1,1], "_", subSample[currAss,2], ".sh", sep="")
    totalStringsSH <- paste("#$ -cwd
#$ -V
#$ -pe smp 1            
#$ -l h_rt=48:0:0      
#$ -l h_vmem=20G

/usr/java/jdk1.7.0_04/bin/java -Xmx9G -jar bin/GenomeAnalysisTK.jar \\
-I", bamFile,"\\
-R /data/BCI-EvoCa/william/referenceHG19/hs37d5.fa \\
-T DepthOfCoverage \\
-o", outFile,"\\
-L", bedNameApo,"\\
--start 1 \\
--stop 1000 \\
--nBins 999
")
    lapply(totalStringsSH, write, scriptName, append=FALSE)
  }  
    
}

