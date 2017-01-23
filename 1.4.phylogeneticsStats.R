# script reads in list of tree files and outputs table of stats

################  notes  ##################
#
# 1. get most parsimonious tree from given sets 
#                 |
#                 v
# 2. using Penny and Hendy (1985) distance, make comparison stats table
#                 |
#                 v
# 3. produce nearest taxa comparison stats
#                 |
#                 v
# 4. test leaf length distributions and plot
#
#
##############   funtions   ################

#compare two trees by phylogenetic taxa distance and report number of differences
#tree1 <- totalTreeList[[currComp]][[1]];tree2 <- tempTree; noTaxaSub <- noTaxa; setNamesSub<- setNames
compareTaxa <- function(tree1, tree2, noTaxaSub, setNamesSub){
  #temp comparison table genome tree vs random (generated below)
  tempNearestTaxaTab <- data.frame(matrix(NA, nrow=noTaxaSub, ncol=2))
  names(tempNearestTaxaTab) <- c("tree1", "tree2")
  row.names(tempNearestTaxaTab) <- setNamesSub
  
  tempCophy <- as.list(NA)
  tempCophy[[1]] <- cophenetic(tree1)
  tempCophy[[2]] <- cophenetic(tree2)
  
  #get nearest taxa distances
  for(currCophy in 1:2){
    for(currName in 1:length(setNamesSub)){
      sampleDist <- tempCophy[[currCophy]][, setNamesSub[currName]]
      sampleDist <- sampleDist[sampleDist!=0]
      sampleDistNames <- names(sampleDist)[which(min(sampleDist) == sampleDist)]
      if(length(sampleDistNames)==1){
        tempNearestTaxaTab[setNames[currName], currCophy] <- sampleDistNames 
      }else{
        tempNearestTaxaTab[setNames[currName], currCophy] <- paste(sampleDistNames, collapse = ":")
      }
    }
  }
  taxaComparisonTemp <- table(tempNearestTaxaTab[[1]] == tempNearestTaxaTab[[2]])
  return(taxaComparisonTemp)
}



############## main program ################

library(apTreeshape)
library(ape)
library(phangorn)
library(irr)

#input sampleList from commandline arguments
arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=3){
	stop("\n#### please use syntax > Rscript 2.2.phylogeneticStats.R < sample list file > < holding directory > < prepended.name > ####\n")
}
sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
holdingDir <- arguments[2]
namePrepended <- arguments[3]

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.Set.10.filt.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- "2.phylogenetics/"
namePrepended <- ".tre"
namePrepended2 <- ".allTreeSearch.tre"

SNVprep <- ".snv.annoVar.variant_function.txt"
#SNVprep <- ".snv.annoVar.exonic_variant_function.0.01.txt"
SNVdir <- "1.platypusCalls/somaticTotal/"

#get sample names
sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[-c(20:23)]

noSets <- length(sampleNames)

treeTypes <- c("somaticTotal.0.01", "exome.0.01", "nonSyn.0.01", "synon.0.01", "CNVphylo", "somaticTotalCNVfiltered")
comparisonFileName <- "phyloTreeComparison.all.csv"


######## 1.loop through sample sets and search trees to get the most parsimonious and shape stats ########
for(eventType in 1:length(treeTypes)){
  for(i in 1:noSets){
    subList <- subset(sampleList, sampleList[1]== sampleNames[i])
    
    treFile <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], namePrepended, sep="")
    
    treFileTotal <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], namePrepended2, sep="")
    normalIndex <- subList[1,7]+1
    normalName <- subList[normalIndex,2]
    histFile <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], ".hist.txt", sep="")
    
    treeList <- as.list(0)
    
    #get specific tree set from file
    if(file.exists(treFileTotal)){
      treeList <- read.nexus(file=treFileTotal)
      histTab <- read.table(file=histFile, header=FALSE, sep="\t")
    }else if(file.exists(treFile)){
      treeList <- read.nexus(file=treFile)
      histTab <- NA
    }else{
      next
    }
    
    #no of trees
    notrees <- length(treeList)
    
    #setup results table
    resultsTable <- as.data.frame(matrix(NA, nrow=notrees, ncol=3))
    names(resultsTable) <- c("tree length", "colless test (yule)", "colless test (PDA)")
    row.names(resultsTable) <- paste("tree", c(1:notrees))
    
    #for each tree get stats
    for(j in 1:notrees){
      print(paste("calculating stats for tree", j, "of", treeTypes[eventType], subList[1,1])) 
      
      #convert to single tree object
      phyloTree <- treeList[[j]]
      
      #resolve polytomies
      phyloTree <-multi2di(phyloTree, random = TRUE)
      
      #get tree length
      resultsTable[j,1] <- sum(phyloTree$edge.length)
      
      #drop normal sample
      treeTemp <- drop.tip(phyloTree, normalName)
      
      #convert to treeshape object
      phyloShape <- as.treeshape(treeTemp)
      
      #perform stats (silently)
      dummyFile <- file()
      sink(file=dummyFile)
      if(length(phyloShape$names) < 5){
        yuleTemp <- colless.test(phyloShape, model="yule", n.mc=1000, alternative="greater")
      }else{
        yuleTemp <- likelihood.test(phyloShape, model="yule", alternative = "greater")
      }
      pdaTemp <- colless.test(phyloShape, model="pda", n.mc=1000, alternative="greater")
      sink()
      close(dummyFile)
      
      #populate table
      resultsTable[j,2] <- yuleTemp$p.value
      resultsTable[j,3] <- pdaTemp$p.value
    }
    
    #reorder table
    resultsTable <- resultsTable[order(resultsTable[1]), ]
    #parseTreeList[i, treeTypes[eventType]] <- row.names(resultsTable)[1]
    
    write.table(resultsTable, file=paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1], "/", subList[1,1], ".treeStats.txt", sep=""), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    
    #output most parsimonous tree
    parseFile <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], ".parseTree.tre", sep="")
    mostParseTree <- as.numeric(strsplit(rownames(resultsTable)[1], " ")[[1]][2])
    parseTree <- treeList[[mostParseTree]]
    write.tree(parseTree, file=parseFile)
    
    #plot tree distribution file
    if(length(histTab)>1){
      pdf(file=paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], ".tree_dist.pdf", sep=""), onefile=TRUE, width=5, height=5)
      barplot(histTab[[2]], main="tree length distribution", xlab="tree length", ylab="frequency", xaxt="n", space=0)
      lengthInt <- c(min(histTab[1]):max(histTab[1]))
      axis(side=1, at=c(0.5:(length(histTab[[1]])-0.5)), labels=lengthInt, las=2, tick=FALSE)
      dev.off()
    }
    
    #plot graphs of distribution
    pdf(file=paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], ".tree_p.value.pdf", sep=""), onefile=TRUE, width=5, height=5)
      plot(resultsTable[[2]], resultsTable[[1]], main="tree length vs p-value from yule", xlab="p-value from yule", ylab="tree length")
    dev.off()
  }
}



######## 2.loop through tree types and populate stats table using Penny and Hendy (1985) distance method ########
statsPHtab <- data.frame(matrix(NA, ncol=6, nrow=length(sampleNames)))
names(statsPHtab) <- c("exome", "nonSyn", "synon", "CNVphylo", "diploid", "nTaxa")
row.names(statsPHtab) <- sampleNames

statsBranchSum <- data.frame(matrix(NA, ncol=6, nrow=length(sampleNames)))
names(statsBranchSum) <- c("exome", "nonSyn", "synon", "CNVphylo", "diploid", "nTaxa")
row.names(statsBranchSum) <- sampleNames

totalTreeList <- as.list(NA)
for(currSamp in 1:length(sampleNames)){
  subList <- subset(sampleList, sampleList[1]==sampleNames[currSamp])
  setNames <- subList[[2]]
  normName <- subList[subList[1,7] + 1, 2]
  setNames <- setNames[setNames != normName]
  noTaxa <- length(setNames)
  
  #summary results table
  sumTreeComp <- data.frame(matrix(NA, ncol=length(treeTypes), nrow=length(treeTypes)))
  names(sumTreeComp) <- treeTypes
  row.names(sumTreeComp) <- treeTypes
  
  statsBranchDist <- data.frame(matrix(NA, ncol=length(treeTypes), nrow=length(treeTypes)))
  names(statsBranchDist) <- treeTypes
  row.names(statsBranchDist) <- treeTypes
  
  #get each tree into a list
  treeList <- as.list(NA)
  for(currTree in 1:length(treeTypes)){
    nameTemp <- paste(subList[1,6], holdingDir, treeTypes[currTree], "/", subList[1,1],"/", subList[1,1], ".parseTree", namePrepended, sep="")
    if(!file.exists(nameTemp)){
      print(paste("removed phylogenetic set from analysis of", sampleNames[currSamp], treeTypes[currTree]))
      treeList[[currTree]] <- NULL
      next
    }
    phyloTree <- read.tree(nameTemp)
    
    #resolve polytomies
    phyloTree <- multi2di(phyloTree, random = TRUE)
    
    #drop normal sample and save to list
    #treeList[[currTree]] <- drop.tip(phyloTree, normName)
    treeList[[currTree]] <- phyloTree
    names(treeList)[[currTree]] <- treeTypes[currTree]
  }
  totalTreeList[[currSamp]] <- treeList
  names(totalTreeList)[[currSamp]] <- subList[1,1]
  
  #generate all trees for given number of taxa that will be used in p.value generation. Penny and Hendy (1985) distance method used as tree comparison method
  allTreeTopologies <- allTrees(noTaxa, tip.label=setNames)
  saveVal <- c()
  for(currPlot in 1:length(allTreeTopologies)){
    #total distance measure for noTaxa number of trees, 
    testTree <- root(phy = allTreeTopologies[[currPlot]], outgroup = normName, resolve.root = TRUE)
    saveVal[currPlot] <- dist.topo(as.phylo(as.treeshape(treeList[[1]])), as.phylo(as.treeshape(testTree)), method="PH85")
  }
  #remove self comparison
  testTabDist <- table(saveVal)
  testTabDist <- testTabDist[testTabDist != 0]

  #now get actual Penny and Hendy distance for data
  for(currX in 1:length(sumTreeComp)){
    for(currY in 1:length(sumTreeComp)){
      if(currX == currY){
        next
      }else{
        if(is.null(treeList[[currY]]) | is.null(treeList[[currX]])){
          next
        }
        actualDist <- dist.topo(treeList[[currX]], treeList[[currY]], method="PH85")
        statsBranchDist[currX, currY] <- actualDist
        #get pvalue
        if(sum(testTabDist[as.numeric(names(testTabDist)) > actualDist]) == 0){
          #set pvalue to 1 as difference is certain
          pvalueSave <- 1
        }else{
          #get true p.value
          pvalueSave <- sum(testTabDist[as.numeric(names(testTabDist)) <= actualDist]) / sum(testTabDist) 
        }
        sumTreeComp[currX, currY] <- pvalueSave
      }
    }
  }
  #write table of stats
  statTabName <- paste(subList[1,6], holdingDir, "comparisonStats/", subList[1,1], ".", noTaxa, ".comparisonStats.txt", sep="")
  write.table(sumTreeComp, file=statTabName, sep="\t", quote = FALSE, col.names = TRUE, row.names=TRUE)
  
  statsBranchSum[currSamp, c(1:5)] <- statsBranchDist[1, c(2:6)]
  statsBranchSum[currSamp, 6] <- noTaxa
  
  statsPHtab[currSamp, c(1:5)] <- sumTreeComp[1, c(2:6)]
  statsPHtab[currSamp, 6] <- noTaxa
} 

#write Penny and Hardy distance p-values
PHTabName <- paste(subList[1,6], holdingDir, "comparisonStats/P&HcomparisonStats.txt", sep="")
write.table(statsPHtab, file=PHTabName, sep="\t", quote = FALSE, col.names = TRUE, row.names=TRUE)

#write Penny and Hardy distance
PHTabNameDist <- paste(subList[1,6], holdingDir, "comparisonStats/P&Hdistance.txt", sep="")
write.table(statsBranchSum, file=PHTabNameDist, sep="\t", quote = FALSE, col.names = TRUE, row.names=TRUE)




######## 3.produce nearest taxa comparisons ########

#note: p-values represents probability of obtaining a comparison with that many taxa agreements by chance

#make output table
statsTab <- data.frame(matrix(NA, nrow=noSets, ncol=(length(treeTypes)-1)))
names(statsTab) <- treeTypes[2:6]
row.names(statsTab) <- sampleNames

for(currComp in 1:length(sampleNames)){
  print(paste("#### producing stats for sample set", sampleNames[currComp]))
  
  #use totalTreeList produced from above analysis
  subList <- subset(sampleList, sampleList[1]==sampleNames[currComp])
  setNames <- subList[[2]]
  normName <- subList[subList[1,7] + 1, 2]
  setNames <- setNames[setNames != normName]
  noTaxa <- length(setNames)
  
  #make holding tables
  statsTree <- as.list(NA)
  statsDistList <- c()
  
  #generate and compare random trees
  for(currStats in 1:10000){
    #template tree to randomly alter
    statsTree[[currStats]] <- totalTreeList[[currComp]][[1]]
    
    #get random taxa order
    branchLengths <- statsTree[[currStats]]$edge.length
    branchOrder <- sample(c(1:length(branchLengths)), size=length(branchLengths))
    statsTree[[currStats]]$edge.length <- branchLengths[branchOrder]
    
    #get random branch length order
    taxaNames <- statsTree[[currStats]]$tip.label
    taxaOrder <- sample(c(1:length(taxaNames)), size=length(taxaNames))
    statsTree[[currStats]]$tip.label <- taxaNames[taxaOrder]
    
    #compare trees and get number of nearest taxa differences
    taxaCompDist <- compareTaxa(totalTreeList[[currComp]][[1]], statsTree[[currStats]], noTaxa, setNames)
    
    if("TRUE" %in% names(taxaCompDist)){
      statsDistList[currStats] <- as.numeric(taxaCompDist["TRUE"])
    }else{
      statsDistList[currStats] <- 0
    }
  }
  #convert stats to table
  statsDistListTab <- table(statsDistList)
  
  #generate real comparisons across mutation types and compare to random distribution
  counter <- 1
  for(currTree in 2:length(treeTypes)){
    tempTree <- totalTreeList[[currComp]][[currTree]]
    if(is.null(tempTree)){
      next
    }else{
      taxaCompDist <- compareTaxa(totalTreeList[[currComp]][[1]], tempTree, noTaxa, setNames)
    }
    if("TRUE" %in% names(taxaCompDist)){
      currDist <- as.numeric(taxaCompDist["TRUE"])
    }else{
      currDist <- 0
    }
    sumRight <- sum(statsDistListTab)
    if(currDist %in% as.numeric(names(statsDistListTab))){
      pPoint <- which(currDist == as.numeric(names(statsDistListTab)))
      sumLeft <- sum(statsDistListTab[pPoint:length(statsDistListTab)])
      simPval <- sumLeft / sumRight
    }else{
      simPval <- 1/sumRight
    }
      
    #add to results table
    statsTab[currComp, counter] <- simPval
    counter <- counter + 1
  }
}

#write results table
statsTable <- paste(subList[1,6], holdingDir, "comparisonStats/taxaDistance.total.txt", sep="")
write.table(statsTab, file=statsTable, sep="\t", row.names=TRUE, col.names = TRUE)



######## 4. test leaf length distributions and plot ########

#clonal proportion table
clonalTab <- data.frame(matrix(NA, ncol=3, nrow=length(sampleNames)))
rownames(clonalTab) <- sampleNames
names(clonalTab) <- c("clonal", "subclonal", "regional")

#prop vector
leafPropVect <- c()
counter <- 1

#use leaf lengths to test against a normal distribution as evidence for selection
for(currSam in 1:length(sampleNames)){
  
  currentSample <- subset(sampleList, sampleList[1]==sampleNames[currSam])
  normalName <- currentSample[currentSample[1,7]+1, 2]
  currSamples <- currentSample[[2]]
  currSamples <- currSamples[currSamples!=normalName]
  
  SNVfileName <- paste(currentSample[1,6], SNVdir, currentSample[1,1], "/", currentSample[1,1], SNVprep, sep="")
  SNVdata <- read.table(file=SNVfileName, sep="\t", stringsAsFactors = FALSE, header = FALSE)
  SNVdata <- SNVdata[-1]
  names(SNVdata) <- c("type", "gene", "chrom", "pos", "pos2", "ref", "alt", paste(currentSample[[2]], ".NR", sep=""), paste(currentSample[[2]], ".NV", sep=""))
  
  #make row lists
  leafVect <- c()
  leafCounter <- 1
  branchVect <- c()
  branchCounter <- 1
  trunkVect <- c()
  trunkCounter <- 1
  delRow <- c()
  delCounter <- 1
  
  #assign rows to lists, del row for mutations to be removed
  for(currVar in 1:nrow(SNVdata)){
    assessRow <- SNVdata[currVar, paste(currSamples, ".NV", sep="")]
    assessTab <- table(assessRow!=0)
    if("TRUE" %in% names(assessTab)){
      if(as.numeric(assessTab["TRUE"])==1){
        leafVect[leafCounter] <- currVar
        leafCounter <- leafCounter + 1
      }else if(as.numeric(assessTab["TRUE"])==length(currSamples)){
        trunkVect[trunkCounter] <- currVar
        trunkCounter <- trunkCounter + 1
      }else{
        branchVect[branchCounter] <- currVar
        branchCounter <- branchCounter + 1
      }
    }else{
      delRow[delCounter] <- currVar
      delCounter <- delCounter + 1
    }
  }
  
  #leaf variants table
  leafData <- SNVdata[leafVect, ]
  
  #branch variants
  branchData <- SNVdata[branchVect, ]
  
  #clonal variants
  trunkData <- SNVdata[trunkVect, ]
  
  clonalTab[currSam, ] <- c(nrow(trunkData), nrow(branchData), nrow(leafData)) 
  
  #remove blank rows from main data
  if(!is.null(delRow)){
    SNVdata <- SNVdata[-delRow, ] 
  }
  totalVariants <- nrow(SNVdata)
  
  propTab <- data.frame(matrix(NA, nrow=1, ncol=length(currSamples)))
  names(propTab) <- currSamples
  for(currProp in 1:length(currSamples)){
    currName <- currSamples[currProp]
    propTab[1, currProp] <- nrow(leafData[leafData[[paste(currName, ".NV", sep="")]] > 0, ]) / totalVariants
    leafPropVect[counter] <- propTab[1, currProp]
    names(leafPropVect)[counter] <- currName
    counter <- counter + 1
  }
}


#finalize and save clonal prop graph
clonalTab[4] <- NA
names(clonalTab)[4] <- "total"
for(currTotal in 1:nrow(clonalTab)){
  clonalTab[currTotal, 4] <- sum(clonalTab[currTotal, c(1:3)])
  clonalTab[currTotal, 1] <- clonalTab[currTotal, 1] / clonalTab[currTotal, 4]
  clonalTab[currTotal, 2] <- clonalTab[currTotal, 2] / clonalTab[currTotal, 4] 
  clonalTab[currTotal, 3] <- clonalTab[currTotal, 3] / clonalTab[currTotal, 4]
}
clonalFile <- paste(currentSample[1,6], "2.phylogenetics/comparisonStats/clonalitytable.csv", sep="")
write.csv(clonalTab, file=clonalFile, quote = FALSE, row.names = TRUE)

clonalTab <- clonalTab[nrow(clonalTab):1,]

#plot clonality figure
outputLoc <- paste(currentSample[1,6], "2.phylogenetics/comparisonStats/clonalityProportions.pdf", sep="")
pdf(file=outputLoc, width=5, height=5)
par(mar=c(5,8,2,2), xpd = TRUE, cex.main = 2.5) 
  barplot(as.matrix(t(clonalTab[1:3])), col = c("steelblue", "goldenrod", "salmon"), las=2, horiz = TRUE)
dev.off()

#plot leaf proportion histograms - cancer
outputLoc <- paste(currentSample[1,6], "2.phylogenetics/comparisonStats/cancerLeafHist.pdf", sep="")
pdf(file=outputLoc, width=5, height=5)
par(mar=c(5,5,5,5), xpd = TRUE) 
  shapStat <- round(shapiro.test(leafPropVect[1:73])$p.value, digits = 6) # null = from a normal distrubution 
  hist(leafPropVect[1:73], breaks = seq(0,0.2,0.01), main="cancer leaf proportion distributions", col="salmon")
  text(y=10, x=0.2, labels = shapStat)
dev.off()


#plot leaf proportion histograms - adenoma
outputLoc <- paste(currentSample[1,6], "2.phylogenetics/comparisonStats/adenomaLeafHist.pdf", sep="")
pdf(file=outputLoc, width=5, height=5)
par(mar=c(5,5,5,5), xpd = TRUE) 
  shapStat <- round(shapiro.test(leafPropVect[74:length(leafPropVect)])$p.value, digits = 6) # null = from a normal distrubution 
  hist(leafPropVect[74:length(leafPropVect)], breaks = seq(0,0.2,0.01), main="adenoma leaf proportion distributions", col="salmon")
  text(y=3, x=0.2, labels = shapStat)
dev.off()

