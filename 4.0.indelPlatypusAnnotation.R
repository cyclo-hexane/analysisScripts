# takes indels called from platypus and annotates using intermediate file to overcome annoVar problem with calling deletions

################### notes ###################
#
# call variants using Platypus and Scalpel (perl scripts)
#
# 4.0.indelPlatypusAnnotation (correct columns and annotate using annovar)
#       |
#       v
# 4.1.indelCompileScalpel (compileTables for Scalpel calls)
#       |
#       v
# 4.2.indelConformance (compare Platypus to Scalpel and get conformance table) 
#       |
#       v
# 4.3.phylogeneticsPrep.indels (make trees using PAUP)
#
#
################# libraries #################



################# subroutines #################

makeAnnovar <- function(input, dirLocSub){
  #run annovar program on new file
  system(command=paste("~/bin/annovar/annotate_variation.pl -out " ,dirLocSub ," -build hg19 ", input, " ~/bin/annovar/humandb/",sep=""))
  
  #annovar output files
  annoExoFunct <- paste(dirLocSub,".exonic_variant_function", sep="")
  annoExoFunctNew <- paste(dirLocSub,".exonic_variant_function.txt", sep="")
  
  annoVarFunct <- paste(dirLocSub,".variant_function", sep="")
  annoVarFunctNew <- paste(dirLocSub,".variant_function.txt", sep="")
  
  #rename files with .txt prepend
  system(command=paste("mv ", annoExoFunct, " ", annoExoFunctNew, sep=""))
  system(command=paste("mv ", annoVarFunct, " ", annoVarFunctNew, sep=""))
}




############### main program ################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.indels.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[[1]])

#indel names
holdingDir <- "1.platypusCalls/indels/"

indelName <- ".indel.somatic.txt"

for(currSam in 1:length(sampleNames)){
  print(paste("#### filtering sample ", sampleNames[currSam], " ####",sep=""))
  
  subSample <- subset(sampleList, sampleList[1]==sampleNames[currSam])
  
  #remove normal
  normalSample <- subSample[subSample[1,7]+1, 2]
  
  #get current indel vcf (.txt) file
  indelFile <- paste(subSample[1,6], holdingDir, sampleNames[currSam], "/", sampleNames[currSam], indelName, sep="")
  indelDataIn <- read.table(file=indelFile, header=FALSE, stringsAsFactors=FALSE, sep="\t")
  
  #insert type column
  indelDataIn <- data.frame(append(x=indelDataIn, values=NA, after=1)) 
  
  #insert second position column
  indelDataIn <- data.frame(append(x=indelDataIn, values=NA, after=5))
  indelDataIn[[2]] <- indelDataIn[[3]]
  
  #convert to annotation file
  for(currConv in 1:nrow(indelDataIn)){
    refStrings <- as.character(indelDataIn[currConv, 4])
    altStrings <- as.character(indelDataIn[currConv, 5])
    if(nchar(refStrings) > nchar(altStrings)){
      indelDataIn[currConv, 6] <- "del"
      if(nchar(refStrings) > 2){
        indelDataIn[currConv, 3] <- indelDataIn[currConv, 3] + (nchar(refStrings)-1)
      }
    }else{
      indelDataIn[currConv, 6] <- "ins"
      if(nchar(refStrings) > 2){
        indelDataIn[currConv, 3] <- indelDataIn[currConv, 3] + (nchar(refStrings)-1)
      }
    }
  }
  
  #output new table
  indelFile <- paste(subSample[1,6], holdingDir, sampleNames[currSam], "/", sampleNames[currSam], ".indel.temp.txt", sep="")
  write.table(indelDataIn, file=indelFile, sep="\t", quote = FALSE, col.names=FALSE, row.names=FALSE)
  
  dirLoc <- paste(subSample[1,6], holdingDir, sampleNames[currSam],"/", sampleNames[currSam], sep="")
  makeAnnovar(indelFile, dirLoc)
  
  #write column names into table
  newAnnoIn <- paste(dirLoc,".exonic_variant_function.txt", sep="")
  tabIn <- read.table(file=newAnnoIn, sep="\t", header=FALSE, stringsAsFactors = FALSE)
  tabIn <- tabIn[-1]
  names(tabIn) <- c("type", "gene", "chrom", "posStart", "posEnd", "ref", "alt", "indel", paste(subSample[[2]], ".NR", sep=""), paste(subSample[[2]], ".NV", sep=""))
  
  #remove normals
  remCol1 <- which(paste(normalSample, ".NV", sep="") ==names(tabIn))
  remCol2 <- which(paste(normalSample, ".NR", sep="") ==names(tabIn))
  tabIn <- tabIn[-c(remCol1, remCol2)]
  
  for(currRow in 1:nrow(tabIn)){
    geneTemp <- strsplit(tabIn[currRow, "gene"], ":")
    tabIn[currRow, "gene"] <- paste(geneTemp[[1]][1], ":", geneTemp[[1]][3], ":", geneTemp[[1]][4], sep="")
  }
  
  #write new table
  write.table(tabIn, file=newAnnoIn, sep="\t", quote = FALSE, col.names=TRUE, row.names=FALSE)
}