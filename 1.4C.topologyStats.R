# script reads in list of tree files and outputs table of stats

################  notes  ##################
#
# compare topology of:
# genome - exome
# genome - nonSyn
# genome - synon
#
# output results on table
#
##############   funtions   ################

############## main program ################

library(apTreeshape)
library(ape)
library(phangorn)
library(irr)

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

holdingDir <- "2.phylogenetics/"
namePrepended <- ".allTreeSearch.tre"
namePrepended2 <- ".parseTree.tre"

#get sample names
sampleNames <- unique(sampleList[[1]])
sampleNames <- sampleNames[c(1:11,24:26)]
#sampleNames <- sampleNames[c(12:16)]


noSets <- length(sampleNames)

#make output table
treeTypes <- c("exome.0.01", "nonSyn.0.01", "synon.0.01")
statsTab <- data.frame(matrix(NA, nrow=length(sampleNames), ncol=(length(treeTypes)+1)))
rownames(statsTab) <- sampleNames
names(statsTab) <- c(treeTypes, "nTaxa")

#tree list files
treeLists <- ".treeStats.txt"

######## 1.loop through sample sets and search trees to get the most parsimonious and shape stats ########
pdf(file="~/PhD/CRCproject/2.phylogenetics/comparisonStats/treeTopologies.adenomas.pdf", width = 10, height = (noSets*2.5))
par(mfrow=c(noSets, 4), mar=c(5,5,5,5), xpd=TRUE)
for(i in 1:noSets){
  subList <- subset(sampleList, sampleList[1]== sampleNames[i])
  normName <- subList[(subList[1,7]+1), 2]
  
  #get somaticTotal list and most parsimoious comaprison tree
  somTreeList <- paste(subList[1,6], holdingDir, "somaticTotal.0.01/", subList[1,1],"/", subList[1,1], treeLists, sep="")
  somaticTotalList <- read.table(file=somTreeList, sep="\t", stringsAsFactors = FALSE)
  parsimTree <- as.numeric(strsplit(rownames(somaticTotalList)[1], split = " ")[[1]][2])
  comparisonTreeFile <- paste(subList[1,6], holdingDir, "somaticTotal.0.01/", subList[1,1],"/", subList[1,1], namePrepended, sep="")
  if(!file.exists(comparisonTreeFile)){
    next()
  }
  compTree <- read.nexus(file=comparisonTreeFile)
  compTree <- compTree[[parsimTree]]
  #compTree <- drop.tip(compTree, normName)
  #compTree <- unroot(compTree)
  compTree <- as.treeshape(compTree, model = "pda")
  plot(compTree, main=paste(subList[1,1], "totalData (comparison)"))
  
  
  #loop through each event type
  for(eventType in 1:length(treeTypes)){
    
    #get most parsimonious trees for this event type
    treeStatsFile <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], treeLists, sep="")
    if(file.exists(treeStatsFile)){
      treeStats <- read.table(file=treeStatsFile, sep="\t", stringsAsFactors = FALSE)
      treeStats <- treeStats[treeStats[[1]]==min(treeStats[["tree.length"]]), ]
    }else{
      next()
    }
    
    treFileTotal <- paste(subList[1,6], holdingDir, treeTypes[eventType], "/", subList[1,1],"/", subList[1,1], namePrepended, sep="")
    if(!file.exists(treFileTotal)){
      #exlude from analysis
      next()
    }else{
      eventTypeTree <- read.nexus(file=treFileTotal)
    }
    
    #if there is only one parsimonious tree compare 1:1
    if(nrow(treeStats) == 1){
      treeindex <- as.numeric(strsplit(rownames(treeStats), split = " ")[[1]][2])
      tempTree <- eventTypeTree[[treeindex]]
      #eventTypeTemp <- drop.tip(eventTypeTemp, normName)
      #eventTypeTemp <- unroot(eventTypeTemp)
      tempTree <- as.treeshape(tempTree, model="pda")
      plot(tempTree, main=paste(subList[1,1], treeTypes[eventType]))
      
      #compare tree topologies and record result
      if(all.equal(compTree, tempTree, names = TRUE, use.edge.length = FALSE)){
        statsTab[i, eventType] <- "1/1"
      }else{
        statsTab[i, eventType] <- "0/1"
      }
    }else{
      #compare each tree and record results
      noTrees <- nrow(treeStats)
      noMatchTrees <- 0
      printFlag <- 0
      for(currTree in 1:nrow(treeStats)){
        treeindex <- as.numeric(strsplit(rownames(treeStats), split = " ")[[currTree]][2])
        tempTree <- eventTypeTree[[treeindex]]
        #tempTree <- drop.tip(tempTree, normName)
        #tempTree <- unroot(tempTree)
        tempTree <- as.treeshape(tempTree, model = "pda")
        if(all.equal(compTree, tempTree, names = TRUE, use.edge.length = FALSE)){
          noMatchTrees <- noMatchTrees + 1
          printFlag <- 1
          plot(tempTree, main=paste(subList[1,1], treeTypes[eventType]))
        }
        
        if(printFlag == 0 & currTree==nrow(treeStats)){
          #plot first tree to compare
          tempPlotHolder <- eventTypeTree[[1]]
          #tempPlotHolder <- drop.tip(tempPlotHolder, normName)
          #tempPlotHolder <- unroot(tempPlotHolder)
          tempPlotHolder <- as.treeshape(tempPlotHolder, model = "pda")
          plot(tempPlotHolder, main=paste(subList[1,1], treeTypes[eventType]))
        }
      }
      statsTab[i, eventType] <- paste(noMatchTrees, "/", noTrees, sep="")
      #legend(x=0, y=0, legend = c(statsTab[i, eventType]), cex=0.6)
    }
  }
}
dev.off()

for(currAdd in 1:nrow(statsTab)){
  subList <- subset(sampleList, sampleList[1]== sampleNames[currAdd])
  statsTab[currAdd, "nTaxa"] <- nrow(subList)-1
}

#save table
write.table(statsTab, file="~/PhD/CRCproject/2.phylogenetics/comparisonStats/treeExactMatchStats.adenoma.txt", sep="\t", quote=FALSE, row.names = TRUE, col.names = TRUE)


