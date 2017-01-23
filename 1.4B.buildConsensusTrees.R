# script reads in list of tree files and outputs consensus tree

################  notes  ##################
#
#     1.1.phylogeneticsPrep.annoVar.R
#                 |
#                 V
#     1.4.phylogeneticsStats.R
#                 |
#                 V
#     1.4B.buildConsensusTrees.R
#
#
##############   funtions   ################


############## main program ################

library(apTreeshape)
library(ape)
library(phangorn)
library(irr)

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- "2.phylogenetics/"
namePrepended <- ".tre"
treeStats <- ".treeStats.txt"

#get sample names
sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[-15]

noSets <- "Set.06"

######## 1.loop through sample sets and search trees to get the most parsimonious and shape stats ########

for(i in 1:length(noSets)){
  subList <- subset(sampleList, sampleList[1]== noSets[i])
  
  treFile <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], namePrepended, sep="")
  statsFile <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], treeStats, sep="")
  outputLoc <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], ".consensusTree.pdf", sep="")
  
  normalIndex <- subList[1,7]+1
  normalName <- subList[normalIndex,2]

  #get tree set from file
  treeList <- read.nexus(file=treFile)
  
  #trees required for consensus building
  statInfo <- read.table(file=statsFile, header=TRUE, sep="\t")
  
  #minimum tree evo length
  minTreeLength <- min(as.numeric(statInfo[[1]]))
  
  #get min trees to merge
  statInfo <- rownames(statInfo[statInfo[[1]]==minTreeLength, ])
  treeNoVect <- c()
  for(currStat in 1:length(statInfo)){
    treeNoVect[currStat] <- as.numeric(strsplit(statInfo[currStat], " ")[[1]][2])
  }
  
  #get tree list to merge
  consenList <- as.list(NA)
  for(currMerge in 1:length(treeNoVect)){
    consenList[[currMerge]] <- treeList[[treeNoVect[currMerge]]]
  }
  
  #make consensus tree
  pdf(file=outputLoc, width=5, height=5)
    plot(consensus(consenList))
  dev.off()
}



