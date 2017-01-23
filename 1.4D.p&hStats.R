# script reads in list of tree files and outputs table of stats

################  notes  ##################
#
# using Penny and Hendy (1985) distance, make conparison stats table
#                 
#
##############   funtions   ################


############## main program ################

library(apTreeshape)
library(ape)
library(phangorn)
library(irr)

# #input sampleList from commandline arguments
# arguments <- commandArgs(trailingOnly = TRUE)
# if(length(arguments)!=3){
# 	stop("\n#### please use syntax > Rscript 2.2.phylogeneticStats.R < sample list file > < holding directory > < prepended.name > ####\n")
# }
# sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
# holdingDir <- arguments[2]
# namePrepended <- arguments[3]

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterLynchList.filt.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- "2.phylogenetics/"
namePrepended <- ".tre"
namePrepended2 <- ".allTreeSearch.tre"

#get sample names
sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[c(1:11,24:26)]
#sampleNames <- sampleNames[c(12:16)]

noSets <- length(sampleNames)

treeTypes <- c("exome.0.01", "nonSyn.0.01", "synon.0.01")

######## loop through tree types and populate stats table using Penny and Hendy (1985) distance method ########
statsBranchSum <- data.frame(matrix(NA, ncol=4, nrow=length(sampleNames)))
names(statsBranchSum) <- c("exome", "nonSyn", "synon", "nTaxa")
row.names(statsBranchSum) <- sampleNames

#tree list files
treeListName <- ".treeStats.txt"

for(currSamp in 1:length(sampleNames)){
  subList <- subset(sampleList, sampleList[1]==sampleNames[currSamp])
  setNames <- subList[[2]]
  normName <- subList[subList[1,7] + 1, 2]
  setNames <- setNames[setNames != normName]
  noTaxa <- length(setNames)
  
  print(paste("#### sample ", sampleNames[currSamp], " ####",sep=""))
  
  #get somatic total tree for comparison, and drop normal
  somaName <- paste(subList[1,6], holdingDir, "somaticTotal.0.01/", subList[1,1],"/", subList[1,1], ".parseTree.tre", sep="")
  comparisonTree <- read.tree(file = somaName)
  comparisonTree <- drop.tip(comparisonTree, normName)
  comparisonTree <- unroot(comparisonTree)
  
  statsBranchDist <- data.frame(matrix(NA, ncol=length(treeTypes), nrow=2))
  names(statsBranchDist) <- treeTypes
  rownames(statsBranchDist) <- c("differentBranches", "pvalue")

  #generate all trees for given number of taxa that will be used in p.value generation. Penny and Hendy (1985) distance method used as tree comparison method
  if(noTaxa > 10){
    allTreeTopologies <- as.list(NA)
    #ntaxa=10 would produce 2027025 tree shapes
    for(curGet in 1:200000){
      allTreeTopologies[[curGet]] <- rtree(noTaxa, rooted = FALSE, tip.label = setNames, br = NULL)
    }
  }else{
    allTreeTopologies <- allTrees(noTaxa, tip.label=setNames, rooted = FALSE)
  }
  saveVal <- c()
  
  if(noTaxa <= 8){
    noLoop <- 1:length(allTreeTopologies)
  }else{
    noLoop <- sample(x = c(1:length(allTreeTopologies)), size = 100000, replace = FALSE)
  }
  
  #get branch difference measures (PH85)
  for(currPlot in noLoop){
    #total distance measure for noTaxa number of trees, 
    #testTree <- root(phy = allTreeTopologies[[currPlot]], outgroup = normName, resolve.root = TRUE)
    #testTree <- drop.tip(as.phylo(testTree), normName)
    testTree <- allTreeTopologies[[currPlot]]
    #saveVal[currPlot] <- dist.topo(as.phylo(as.treeshape(comparisonTree)), as.phylo(as.treeshape(testTree)), method="PH85")
    saveVal[currPlot] <- dist.topo(comparisonTree, testTree, method = "PH85")
  }
  #remove self comparison
  testTabDist <- table(saveVal)
  testTabDist <- testTabDist[testTabDist != 0]
  
  #get each tree into a list
  newTreeList <- as.list(NA)
  for(currTree in 1:length(treeTypes)){
    
    #get stats file and determine whether tree is unique
    treeStatsFile <- paste(subList[1,6], holdingDir, treeTypes[currTree], "/", subList[1,1],"/", subList[1,1], treeListName, sep="")
    if(file.exists(treeStatsFile)){
      treeStats <- read.table(file=treeStatsFile, sep="\t", stringsAsFactors = FALSE)
      treeStats <- treeStats[treeStats[[1]]==min(treeStats[["tree.length"]]), ]
    }else{
      next()
    }
    
    #get all trees from brute force search
    nameTemp <- paste(subList[1,6], holdingDir, treeTypes[currTree], "/", subList[1,1],"/", subList[1,1], ".allTreeSearch", namePrepended, sep="")
    nameTemp2 <- paste(subList[1,6], holdingDir, treeTypes[currTree], "/", subList[1,1],"/", subList[1,1], namePrepended, sep="")
    if(file.exists(nameTemp)){
      print(paste("getting tree set for", sampleNames[currSamp], treeTypes[currTree]))
      treeList <- read.nexus(file=nameTemp)
    }else{
      treeList <- read.nexus(file=nameTemp2)
    }
    
    #get most parsimonious tree indexes
    newTreeList <- as.list(NA)
    for(currIndex in 1:nrow(treeStats)){
      newTreeList[[currIndex]] <- treeList[[as.numeric(strsplit(rownames(treeStats)[currIndex], split = " ")[[1]][2])]]
    }
    
    
    #test each tree
    branchLengthList <- c()
    pvalueList <- c()
    for(currComp in 1:nrow(treeStats)){
      currTestTree <- newTreeList[[currComp]]
      
      #resolve polytomies and drop normal
      phyloTree <- multi2di(currTestTree, random = TRUE)
      phyloTree <- drop.tip(phyloTree, normName)
      phyloTree <- unroot(phyloTree)
      
      #get P&H distance measure
      #actualDist <- dist.topo(as.phylo(as.treeshape(phyloTree)), as.phylo(as.treeshape(comparisonTree)), method="PH85")
      actualDist <- dist.topo(phyloTree, comparisonTree, method="PH85")
      branchLengthList[currComp] <- actualDist
      
      #get probability of observing such few number of different branches, significant is p<0.05
      if(sum(testTabDist[as.numeric(names(testTabDist)) > actualDist]) == 0){
        #set pvalue to 1 as difference is certain
        pvalueSave <- 1
      }else{
        #get true p.value, significant pvalue indicates a low likelihood of obtaining a tree that has this few branch partitions (i.e similarity)   
        pvalueSave <- sum(testTabDist[as.numeric(names(testTabDist)) <= actualDist]) / sum(testTabDist) 
      }
      
      pvalueList[currComp] <- pvalueSave
    }
    
    statsBranchDist[1, currTree] <- paste(branchLengthList, collapse = ",")
    statsBranchDist[2, currTree] <- paste(pvalueList, collapse = ",")
  }
  
  #add to summary table
  statsBranchSum[currSamp, c(1:4)] <- c(as.character(statsBranchDist[2, ]), noTaxa)
}


#write Penny and Hardy distance
PHTabNameDist <- paste(subList[1,6], holdingDir, "comparisonStats/P&Hdistance.adenoma.txt", sep="")
write.table(statsBranchSum, file=PHTabNameDist, sep="\t", quote = FALSE, col.names = TRUE, row.names=TRUE)


