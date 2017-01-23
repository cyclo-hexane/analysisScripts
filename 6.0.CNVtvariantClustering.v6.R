# script attempts to order SNV and CNV events in a given cancer and output a progression model
# reguires:
# 1. sampleList file (standard columns)
# 2. CNV segmentation file and
# 3. joint SNV calls from platypus
# 4. Driver gene list

################## notes ##################
#   This version removes all SNVs that are not truncal from clustering
#
#   6.0.CNVtvaraintClustering.R (est alpha/beta and make input table: X.matlab.input.txt)
#             |
#             V
#   AnalyseCNVData.m (matlab script on apocrita, output: X.matlab.output.txt)
#             |
#             V
#   6.1.estimateBifurcations.R (used for plot)
#             |
#             V
#   6.2.assessCNVtiming.heatmap.R (make final figure of timings)
#
#
# this script:
#
#           A. filter valid CNV regions 
#                   |
#                   V
#           B. get SNV for each region
#                   |
#                   V
#           C. calculate alpha and beta
#                   |
#                   V
#           D. create model input files (and split by each biopsy)
#
#
################## assumptions ##################


################## subroutines ##################

#### A. filter valid CNV events and assign phylogenetic location ####
#segDataSub <- segCNV;noSamSub <- noSamples
filterCNVs <- function(segDataSub, noSamSub){
  returnList <- as.list(NA)
  
  #assess segmentations if normal, remove
  removeSegs <- c()
  counterSub <- 1
  segDataSub <- as.data.frame(append(segDataSub, x = NA))
  names(segDataSub)[1] <- "phyloLoc"
  for(currSeg in 1:nrow(segDataSub)){
    assessSeg <- table(segDataSub[currSeg, 6:(5+(noSamSub*2))] == 1)
    if("TRUE" %in% names(assessSeg)){
      if(as.integer(assessSeg["TRUE"])==(2*noSamSub)){
        #don't keep, its diploid
        removeSegs[counterSub] <- currSeg
        counterSub <- counterSub + 1
      }else{
        #keep
      }
    }
  }
  #subset to remove normal regions
  if(length(removeSegs)>0){
    segDataSub <- segDataSub[-removeSegs,]
  }
  
  #assess phylogenetic location
  for(currSeg in 1:nrow(segDataSub)){
    #get current CNV seg into table
    assessVector <- segDataSub[currSeg, 6:(5+(noSamSub*2))]
    assessTable <- data.frame(matrix(NA, ncol=2, nrow=noSamSub, byrow = TRUE))
    seqAp <- seq(1, (noSamSub*2), 2)
    for(currAp in 1:nrow(assessTable)){
      assessTable[currAp, 1] <- assessVector[seqAp[currAp]]
      assessTable[currAp, 2] <- assessVector[seqAp[currAp]+1]
    }
    assessTable[3] <- paste(assessTable[[1]], ":", assessTable[[2]], sep="")
    names(assessTable) <- c("minor", "major", "merged")
    
    #assess phylo by uniqueness of states
    stateTest <- table(assessTable[["merged"]])
    noStates <- length(unique(stateTest))
    
    #assess number of states
    if(noStates == 1){
      #trunkal
      segDataSub[currSeg, 1] <- "T"
    }else if("1:1" %in% names(stateTest)){
      if(stateTest["1:1"] == (noSamSub-1)){
        #true leaf event
        segDataSub[currSeg, 1] <- "L"
      }else if(noStates == 2){
        #its a true branch event
        segDataSub[currSeg, 1] <- "B"
      }else{
        #its a bifuricated branch (two branch events in same region)
        segDataSub[currSeg, 1] <- "BiB"
      }
    }else{
      #if no samples have WT state they are ether branch or trunk
      if(noStates == 2){
        #its a bifuricated trunk (one branch events in same region as a trunk)
        segDataSub[currSeg, 1] <- "BiT"
      }else{
        #phylogenetic loc cannot be determined
        segDataSub[currSeg, 1] <- "unknown >2 bifurcations"
      }
    }
  }
  returnList[[1]] <- segDataSub
  returnList[[2]] <- removeSegs
  return(returnList)
}


#### C. get SNV for each region ####
getSNVs <- function(currPosStartSub, currPosEndSub, currChromSub, currSamSub, vcfDataSub){
  #subset by chromosome location
  vcfDataSub <- subset(vcfDataSub, vcfDataSub$chrom==currChromSub)
  keepListSub <- c()
  counterSub <- 1
  
  #get current sample AF col name
  currSampCol <- paste(currSamSub, ".AF", sep="")
  
  #loop though variants and keep those in CNV
  for(currRow in 1:nrow(vcfDataSub)){
    if(vcfDataSub[currRow, "pos"] >= currPosStartSub & vcfDataSub[currRow, "pos"] <= currPosEndSub & vcfDataSub[currRow, currSampCol] > 0){
      #keep variant as it is in CNV region
      keepListSub[counterSub] <- currRow
      counterSub <- counterSub + 1
    }
  }
  vcfDataSub <- vcfDataSub[keepListSub, c("chrom", "pos", "ref", "alt", currSampCol)]
  return(vcfDataSub)
}


#### C. calculate alpha and beta parameters
#tempSNVtableSub <- tempSNVtable;currentStateSub <- currentState; graphOutDirSub <- graphOutDir; celluEstSub <- cellularityEst; covSub <- covData
calcParameters <- function(tempSNVtableSub, currentStateSub, graphOutDirSub, celluEstSub, covSub){
  #get CNV type (trisomy, BLOH etc)
  if(as.numeric(currentStateSub[1]) == 1 & as.numeric(currentStateSub[2]) == 2){
    CNVstateSub <- "trisomy"
    estBetaCenter <- 0.33
    estAlphaCenter <- 0.67
  }else if(as.numeric(currentStateSub[1]) == 2 & as.numeric(currentStateSub[2]) == 2){
    CNVstateSub <- "tetrasomy"
    estBetaCenter <- 0.25
    estAlphaCenter <- 0.5
  }else if(currentStateSub[1] == 0 & currentStateSub[2] == 1){
    CNVstateSub <- "loss"
    estBetaCenter <- 0.5
    estAlphaCenter <- 1
  }else if(currentStateSub[1] == 0 & currentStateSub[2] == 2){
    CNVstateSub <- "BLOH"
    estBetaCenter <- 0.5
    estAlphaCenter <- 1
  }else{
    CNVstateSub <- "polysomy"
    estBetaCenter <- 0.25
    estAlphaCenter <- 0.5
  }
  
  #get alpha and beta (if correct CNV type)
  if(nrow(tempSNVtableSub) >= 50){
    totalVariants <- nrow(tempSNVtableSub)
    estBetaCenter <- estBetaCenter * celluEstSub 
    estAlphaCenter <- estAlphaCenter * celluEstSub 
    
    Afdata <- tempSNVtableSub[[5]]
    
    #cluster allele frequencies using mixture model, repeating if the same mean is found
    clusteredAF <- normalmixEM(Afdata, mu = c(estBetaCenter, estAlphaCenter), arbvar = FALSE, epsilon = 10)
    clustCounter <- 1
    while((abs(clusteredAF$mu[1] - clusteredAF$mu[2]) < 0.1) & clustCounter < 10000){
      clusteredAF <- normalmixEM(Afdata, mu = c(estBetaCenter, estAlphaCenter), arbvar = FALSE, epsilon = 10)
      clustCounter <- clustCounter + 1
    }
    clustFlag <- 1
    if(abs(clusteredAF$mu[1] - clusteredAF$mu[2]) < 0.1){
      #cluster ill-fitted mark stats as NA
      clustFlag <- 0
    }else{
      #get statistics for mixture model
      clusterStats <- clusteredAF$loglik
    }
    
    #get cluster centre means and assign cluster charatcers to alpha or beta
    alpha <- round(totalVariants * (clusteredAF$lambda[2]), digits = 0)
    beta <- round(totalVariants * (clusteredAF$lambda[1]), digits = 0)
    
    alphaMean <- clusteredAF$mu[2]
    betaMean <- clusteredAF$mu[1]
    
    distBeta <- abs(x = betaMean - estBetaCenter) 
    distAlpha <- abs(x = alphaMean - estAlphaCenter) 
    
    #ci defined as 1sd
    ciAlpha <- clusteredAF$sigma[1]
    ciBeta <- clusteredAF$sigma[2]
    ciHighAlpha <- alphaMean + ciAlpha
    ciLowAlpha <- alphaMean - ciAlpha
    ciHighBeta <- betaMean + ciBeta
    ciLowBeta <- betaMean - ciBeta
    
    #plot graph of clusters
    pdf(file=paste(graphOutDirSub, CNVstateSub, ".pdf", sep=""))
    par(xpd=FALSE)
    
      #plot raw VAFs as density and modelled distribution
      hist(Afdata, xlab="VAF", breaks=seq(-0.1,1,0.01), main="clustering VAF distribution", freq = FALSE, xlim = c(-0.1:1))
      
      alphaDensity <- density(rnorm(n = 1000000, mean=alphaMean, sd=ciAlpha), kernel = "gaussian", from = -0.1, to = 1, cut = 0.01)
      alphaDensityY <- alphaDensity$y * (alpha / sum(alphaDensity$y))
      maxAlphaYdensity <- max(alphaDensityY)
      lines(x=alphaDensity$x, y=alphaDensityY, col="red3", pch=2, lwd=2, add=TRUE)
      
      betaDensity <- density(rnorm(n = 1000000, mean=betaMean, sd=ciBeta), kernel = "gaussian", from = -0.1, to = 1, cut = 0.01)
      betaDensityY <- betaDensity$y * (beta / length(betaDensity$y))
      maxBetaYdensity <- max(betaDensityY)
      lines(x=betaDensity$x, y=betaDensityY, col="blue4", pch=2, lwd=2, add=TRUE)
      
      #plot cluster means
      lines(x=c(alphaMean, alphaMean), y=c(0, maxAlphaYdensity), lty=1, pch=2, lwd=3, col="red3")
      lines(x=c(betaMean, betaMean), y=c(0, maxBetaYdensity), lty=1, pch=2, lwd=3, col="blue4")
      
      #plot estimated cluster centres
      lines(x=c(estBetaCenter, estBetaCenter), y=c(0, 10), lty=2, pch=2, lwd=2, col="blue4")
      lines(x=c(estAlphaCenter, estAlphaCenter), y=c(0, 10), lty=2, pch=2, lwd=2, col="red3")
      
      #add CI for clustered dat
      lines(x = c(ciLowAlpha, ciHighAlpha), y=rep(0, 2), lty=1, pch=2, lwd=2, col="red3") 
      lines(x = c(ciLowBeta, ciHighBeta), y=rep(0, 2), lty=1, pch=2, lwd=2, col="blue4")
      
      #add cluster log-likelihood ratios
      text(x=0.8, y=0.5, paste("cluster LLR:", round(clusterStats, digits = 2)), cex=0.6)
      
      #add cluster stats to graph
      #text(x=alphaMean, y=yMax, paste("mean distances: ", round(distBeta, digits = 6), round(distAlpha, digits = 6)), cex=0.6)
    dev.off()
    
    #check if means are within expected sd boundary
    if((distAlpha < ciAlpha & distBeta < ciBeta) & clustFlag == 1){
      clusterDist <- paste(round(distAlpha, digits = 6), ":", round(distBeta, digits = 6), sep="")
    }else if(clustFlag == 0){
      clusterStats <- "N"
    }else if(distAlpha > ciAlpha){
      clusterStats <- "A"
    }else{
      clusterStats <- "B"
    }
    
    #return parameters
    parameterList <- c(alpha, beta, clusterDist, CNVstateSub, alphaMean, betaMean, clusterStats, paste(graphOutDirSub, CNVstateSub, ".pdf", sep=""))
    
  }else{
    #cannot calculate small regions
    parameterList <- c(nrow(tempSNVtableSub), 0, NA, CNVstateSub, NA, NA, NA, NA)
  }
  return(parameterList)
}



################## libraries ##################

library(fpc)
library(mixtools)

################## main program ##################



#global variables
# 1. sampleList file
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.Set.10.filt.csv", header=FALSE, stringsAsFactors=FALSE)

# 2. CNV segmentation file
CNVPrepended <- ".penalty0.95.baf.gt.txt"

# 3. SNV copy state file
variantCallFile <- ".snv.annoVar.variant_function.0.01.txt"

# 4. purity estimates
purityFile <- ".penalty0.95.purity-cna-baf.txt"
workingDir <- "7.CNVcalls.archive/baf.updated/"
SNVworkingDir <- "1.platypusCalls/somaticTotal.0.01/"
analysisDir <- "6.CNVtiming/inputData.v6/"
purityDir <- "13.purityEstimates/"
setNames <- unique(sampleList[[1]])
setNames <- setNames[c(2,3,5,9,10,11,16)]

# segmentations to exclude
excDir <- "6.CNVtiming/exclusionRegions/"
excludeSegs <- c(8)

#store mean distances
distVect <- data.frame(matrix(NA, nrow=0, ncol=1))
names(distVect) <- "distances"

#main program
for(currSet in 1:length(setNames)){
  
  ####### setup variables and get segmentation table
  currentSamples <- subset(sampleList, sampleList[1]==setNames[currSet])
  
  print(paste("starting analysis for set", currentSamples[1,1]))
  
  #make output directories
  system(command=paste("mkdir ", currentSamples[1,6], analysisDir, currentSamples[1,1], sep=""))
  system(command=paste("mkdir ", currentSamples[1,6], analysisDir, currentSamples[1,1], "/exclude/", sep=""))
  system(command=paste("mkdir ", currentSamples[1,6], analysisDir, currentSamples[1,1], "/include/", sep=""))
  
  #remove normal from list
  normalSample <- which(currentSamples[[2]] %in% currentSamples[currentSamples[1,7]+1, 2])
  currentSamples <- currentSamples[-normalSample, ]
  noSamples <- nrow(currentSamples)
  sampleNames <- currentSamples[[2]]
  
  #get cellularity estimates
  celluEst <- read.table(file=paste(currentSamples[1,6], purityDir, currentSamples[1, 1], purityFile, sep=""), sep=" ")
  
  #get CNV segmentations
  segCNVfile <- paste(currentSamples[1,6], workingDir, currentSamples[1,1], CNVPrepended, sep="")
  segCNV <- read.table(file=segCNVfile, sep="\t", fill=TRUE, stringsAsFactors=FALSE, header=TRUE)
  
  #set ambigous calls to WT
  segCNV[segCNV=="X"] <- 1
  
  #save sampled genome length then filter small regions
  segCNV["nloci"] <- (segCNV[["last.locus"]] - segCNV[["first.locus"]]) / 1000000
  totalGenomeMeasured <- sum(segCNV[[3]]) * 1000000
  
  #ensure states are numeric
  for(currChg in 5:ncol(segCNV)){
    segCNV[currChg] <- as.numeric(segCNV[[currChg]])
  }
  
  
  
  ####### A. filter CNVs for truncal and diploid regions  ####### 
  
  returnedInfo <- filterCNVs(segCNV, noSamples)
  filteredSegCNVs <- returnedInfo[[1]]
  if(length(returnedInfo)>1){
    diploidRegions <- segCNV[returnedInfo[[2]], ]
  }else{
    diploidRegions <- NA
  }
  
  #keep only those CNVs that are clonal (unhash for branch as well)
  nonTimedEvents <- filteredSegCNVs[(filteredSegCNVs[["phyloLoc"]]!="T" & filteredSegCNVs[["phyloLoc"]]!="BiT") | filteredSegCNVs[["nloci"]] < 10, ]
  filteredSegCNVs <- filteredSegCNVs[(filteredSegCNVs[["phyloLoc"]]=="T" | filteredSegCNVs[["phyloLoc"]]=="BiT") & filteredSegCNVs[["nloci"]] >= 10, ]
  
  #output non-timed segmentations
  nonTimedEvents <- nonTimedEvents[order(nonTimedEvents[["chr"]], nonTimedEvents[["first.locus"]]), ]
  outFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".nonTimedEvents.txt", sep="")
  write.table(nonTimedEvents, file=outFileName, sep=",", quote=FALSE, row.names=FALSE)
  
  
  
  ####### B. get vcf and SNVs for this set ########
  
  vcfFileName <- paste(currentSamples[1,6], SNVworkingDir, currentSamples[1,1], "/", currentSamples[1,1], variantCallFile, sep="")
  vcfData <- read.table(file=vcfFileName, sep="\t", stringsAsFactors = FALSE, header=FALSE)
  
  #remove normal sample
  vcfData <- vcfData[-c((8 + normalSample) , (8 + (noSamples + 1) + normalSample))]
  
  #calculate allele frequencies and change X/Y to 23/24
  for(y in 1:noSamples){
    vcfData[(ncol(vcfData)+1)] <- as.vector(vcfData[(8+(noSamples)+y)]) / as.vector(vcfData[8+y])
  }
  names(vcfData) <- c("line", "region", "variant", "chrom", "pos", "pos2", "ref", "alt", paste(sampleNames, ".NR", sep=""), paste(sampleNames, ".NV", sep=""), paste(sampleNames, ".AF", sep=""))
  chromData <- vcfData[["chrom"]]
  chromData[chromData=="X"] <- 23
  chromData[chromData=="Y"] <- 24
  vcfData["chrom"] <- chromData
  
  #remove variants not in CNV regions that are hypermutated regions, set to blank if there are no reagons to remove
  if(excludeSegs == currSet){
    excFile <- paste(currentSamples[1,6], excDir, currentSamples[1,1], ".txt", sep="")
    excTab <- read.table(file=excFile, sep="\t", stringsAsFactors = FALSE, header=FALSE)
  }else{
    excTab <- data.frame(0, nrow=0, ncol=3)
  }
  
  vcfDataUnfilt <- vcfData
  
  #remove snv if not truncal or in exclusion table above
  remRow <- c()
  remCounter <- 1
  for(currSNV in 1:nrow(vcfData)){
    SNVchrom <- vcfData[currSNV, 4]
    SNVpos <- vcfData[currSNV, 5]
    
    #if no CNV region exists for this chromosome, remove
    subCNV <- filteredSegCNVs[filteredSegCNVs[["chr"]]==SNVchrom, ]
    if(nrow(subCNV) == 0){
      remRow[remCounter] <- currSNV
      remCounter <- remCounter + 1
      next()
    }
    
    #if there are CNV regions to filter, test current variant
    subHyper <- excTab[excTab[[1]]==SNVchrom, ] 
    if(nrow(subHyper) != 0){
      #if truncal, see if variant is in filter chromosome regions (subHyper)
      for(currHyp in 1:nrow(subHyper)){
        if(SNVpos >= subHyper[currHyp, 2] & SNVpos <= subHyper[currHyp, 3]){
          remRow[remCounter] <- currSNV
          remCounter <- remCounter + 1
          next()
        }
      }
    }
    
    #variant not truncal, remove
    if(as.numeric(table(vcfData[currSNV, paste(sampleNames, ".AF", sep="")] > 0)["TRUE"]) != noSamples){
      remRow[remCounter] <- currSNV
      remCounter <- remCounter + 1
      next()
    }
    
    #finally, see if SNV in a CNV region
    for(currTest in 1:nrow(subCNV)){
      if(SNVpos >= subCNV[currTest, "first.locus"] & SNVpos <= subCNV[currTest, "last.locus"]){
        break()
      }else if(currTest == nrow(subCNV)){
        #remove as not in any CNV
        remRow[remCounter] <- currSNV
        remCounter <- remCounter + 1
      }
    }
  }
  
  vcfData <- vcfData[-remRow, ]
  
  
  
  
  
  ####### C. loop through segmentations and attempt to cluster ####### 
  
  #setup results table
  CNVtotalList <- data.frame(matrix(NA, nrow=(noSamples * nrow(filteredSegCNVs)), ncol=17))
  names(CNVtotalList) <- c("SetName", "SampleID", "CNV.ID", "Chromosome", "CNVlength", "start", "end", "alpha", "beta", "meanDist", "type", "seClust1", "seClust2", "clusterStats", "phyloLoc", "cellularity", "filter")
  CNVtotalList[1] <- currentSamples[1,1]
  CNVtotalList[2] <- currentSamples[[2]]

  CNVcounter <- 1
  colSeqMaj <- paste(sampleNames, "_Major", sep="")
  colSeqMin <- paste(sampleNames, "_Minor", sep="")
  for(currCNV in 1:nrow(filteredSegCNVs)){
    currChrom <- filteredSegCNVs[currCNV, "chr"]
    currPosStart <- filteredSegCNVs[currCNV, "first.locus"]
    currPosEnd <-  filteredSegCNVs[currCNV, "last.locus"]
    CNVlength <- currPosEnd - currPosStart
    
    #for each sample with CNV present get SNVs, calculate alpha and beta and add to table
    for(currSam in 1:noSamples){
      currMajCol <- colSeqMaj[currSam]
      currMinCol <- colSeqMin[currSam]
      currentState <- filteredSegCNVs[currCNV, c(currMinCol, currMajCol)]
      if(filteredSegCNVs[currCNV, "chr"] == 23){
        intervalID <- paste("X:", (filteredSegCNVs[currCNV, "first.locus"]+1), "-", filteredSegCNVs[currCNV, "last.locus"], sep="")
      }else{
        intervalID <- paste(filteredSegCNVs[currCNV, "chr"], ":", (filteredSegCNVs[currCNV, "first.locus"]+1), "-", filteredSegCNVs[currCNV, "last.locus"], sep="")
      }
      
      currID <- paste(currentSamples[currSam, 2], "-chr", currChrom, "-", currentState[1], ".", currentState[2], "-", CNVlength, sep="")
      
      print(paste("##### processing ID:", currID, "#####") )
      
      #current cellularity value
      cellularityEst <- celluEst[celluEst[[3]]==currentSamples[currSam, 2], 4]
      
      # 1. get SNV for each region
      tempSNVtable <- getSNVs(currPosStart, currPosEnd, currChrom, sampleNames[currSam], vcfData)
      
      #get coverage information from biopsy specific file
      covFileName <- paste('~/PhD/CRCproject/6.CNVtiming/gatk_coverage/', setNames[currSet], '_', currentSamples[currSam, 2], '.sample_interval_summary', sep = '')
      covData <- read.table(covFileName, header = T, fill = T)
      covData <- covData[covData[["Target"]]==intervalID, "average_coverage"]
      
      # 2. calculate alpha and beta
      graphOutDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currID, ".", sep="")
      parameters <- calcParameters(tempSNVtable, currentState, graphOutDir, cellularityEst, covData)
      
      #if parameters not possible NA is placed instead of calculated values for alpha and beta
      CNVtotalList[CNVcounter,] <- c(currentSamples[1,1], sampleNames[currSam], currID, currChrom, CNVlength, currPosStart, currPosEnd, parameters[1:7], filteredSegCNVs[currCNV, "phyloLoc"], cellularityEst, "NA")
      
      #check for filter criteria and mark accordingly, also move graph to correct folder
      filterString <- c()
      filtCounter <- 1
      if(is.na(parameters[3])){
        filterString <- "lowVariants"
      }else if(parameters[3] == "A" | parameters[3] == "B"){
        if(parameters[3] == "A"){
          #mean of alpha cluster out of confidence 
          filterString <- "alphaCluster"
        }else{
          #mean of beta cluster out of confidence 
          filterString <- "betaCluster"
        }
        newDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/exclude/", sep="") 
        system(command = paste("mv", parameters[8], newDir))
      }else{
        if(parameters[4] == "polysomy" | parameters[4] == "loss"){
          filterString[filtCounter] <- "copyState"
          filtCounter <- filtCounter + 1
        }
        if(filtCounter != 1){
          newDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/exclude/", sep="") 
          system(command = paste("mv", parameters[8], newDir))
        }else{
          filterString[1] <- "PASS"
          newDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/include/", sep="")
          system(command = paste("mv", parameters[8], newDir))
        }
      }
      
      #add filter to table
      CNVtotalList[CNVcounter, "filter"] <- paste(filterString, collapse = ":")
      
      #update table
      CNVcounter <- CNVcounter + 1
    }
  }
  
  #order by sample
  CNVtotalList <- CNVtotalList[order(CNVtotalList["SampleID"], as.numeric(CNVtotalList[["Chromosome"]]), as.numeric(CNVtotalList[["start"]]), decreasing = TRUE), ]
  
  #save table
  inputFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".input.parameters.total.txt", sep="")
  write.table(CNVtotalList, file=inputFileName, sep=",", quote=FALSE, row.names=FALSE)
  
  #filter for variants with PASS
  CNVtotalList <- CNVtotalList[CNVtotalList[["filter"]]=="PASS",]
  CNVtotalList["filter"] <- NULL
  inputFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".input.parameters.txt", sep="")
  write.table(CNVtotalList, file=inputFileName, sep=",", quote=FALSE, row.names=FALSE)
  
  #save mean distances to vector
  disTemp <- data.frame(matrix(nrow = nrow(CNVtotalList), ncol=1))
  names(disTemp) <- "distances"
  disTemp[1] <- CNVtotalList[["stats"]]
  distVect <- rbind(distVect, disTemp)
  
  #reformat table for model input
  if(nrow(CNVtotalList) > 0){
    #reformat for model input
    for(currCNV in 1:nrow(CNVtotalList)){
      if(CNVtotalList[currCNV, "type"] == "tetrasomy"){
        CNVtotalList[currCNV, "stats"] <- 2
        CNVtotalList[currCNV, "type"] <- 2
      }else if(CNVtotalList[currCNV, "type"] == "trisomy"){
        CNVtotalList[currCNV, "stats"] <- 1
        CNVtotalList[currCNV, "type"] <- 2
      }else if(CNVtotalList[currCNV, "type"] == "BLOH"){
        CNVtotalList[currCNV, "stats"] <- 0
        CNVtotalList[currCNV, "type"] <- 2
      }else{
        CNVtotalList[currCNV, "stats"] <- 0
        CNVtotalList[currCNV, "type"] <- 0
      }
    }
    CNVtotalList["start"] <- NULL
    CNVtotalList["end"] <- NULL
    CNVtotalList[10:11] <- CNVtotalList[6:7]
    CNVtotalList <- CNVtotalList[-c(1,6:7,12:14)]
    names(CNVtotalList) <- c("SetName", "CNV_ID", "Chromosome", "length", "a2", "b2", "alpha", "beta")
    
    #write timing file
    inputFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".matlab.input.txt", sep="")
    write.table(CNVtotalList, file=inputFileName, sep=",", quote=FALSE, row.names=FALSE)
  }
  
  
  
  #get diploid region counts for model parameters (hash diploid regions to gte proper calculation)
  diploidRegions <- 0
  if(length(diploidRegions)==1){
    diploidCounter <- 0
    diploidMB <- 0
  }else{
    diploidCounter <- 0
    for(currCNV in 1:nrow(diploidRegions)){
      currPosStart <- diploidRegions[currCNV, "first.locus"]
      currPosEnd <-  diploidRegions[currCNV, "last.locus"]
      currChrom <- diploidRegions[currCNV, "chr"]
      
      #temp holding vector
      tempSNVcounter <- c()
      
      #count variants in each sample for this diploid region
      for(currSam in 1:length(sampleNames)){
        # B. get SNV for each region
        tempSNVcounter[currSam] <- nrow(getSNVs(currPosStart, currPosEnd, currChrom, sampleNames[currSam], vcfDataUnfilt))
      }
      
      #update counter
      diploidCounter <- diploidCounter + mean(tempSNVcounter)
    }
    diploidMB <- sum(diploidRegions[["nloci"]] * 1000000)
  }
  
  diploidOutput <- data.frame(matrix(c(diploidMB, round(diploidCounter, digits = 0)), byrow = FALSE, nrow=1, ncol=2))
  names(diploidOutput) <- c("length", "mutations")
  
  dipFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".normal.stats.txt", sep="")
  write.table(diploidOutput, file=dipFileName, sep=",", quote=TRUE, row.names=FALSE)
  
} 
  



####### plot mean distances, for inspection ########
distVect <- unlist(strsplit(distVect[[1]], split = ":"))
distPlot <- data.frame(matrix(NA, nrow=(length(distVect)/2), ncol=2))
counter <- 1
for(i in seq(1, length(distVect), 2)){
  distPlot[counter, 1] <- distVect[i]
  distPlot[counter, 2] <- distVect[i+1]
  counter <- counter + 1
}

hist(as.numeric(distPlot[[2]]), breaks = seq(0,1,0.01))
  



####### D. check segmentation plots manually and move segmentations not passed to include folder ####### 

#loop though and reassess include graphs
for(currSet in 1:length(setNames)){
  currentSamples <- subset(sampleList, sampleList[1]==setNames[currSet])
  
  print(paste("re-assessing set", currentSamples[1,1]))
  
  inputClustFile <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".input.parameters.txt", sep="")
  clustIn <- read.table(file=inputClustFile, sep=",", header = TRUE, stringsAsFactors = FALSE)
  
  #get include event list from graphs in include directory
  includeDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/include/", sep="")
  graphsIncluded <- list.files(path = includeDir, all.files = FALSE)
  
  #parse names
  tempFileHold <- matrix(unlist(strsplit(graphsIncluded, "\\.")), ncol = 4, nrow = length(graphsIncluded), byrow = TRUE)
  CNVsIncluded <- paste(tempFileHold[,1], tempFileHold[,2], sep=".")
  
  for(currCNV in 1:length(CNVsIncluded)){
    clustIn[clustIn[["CNV.ID"]]==CNVsIncluded[currCNV], "filter"] <- "PASS"
  }
  
  #filter for variants with PASS
  clustIn <- clustIn[clustIn[["filter"]]=="PASS",]
  clustIn["filter"] <- NULL
  inputFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".input.parameters.AutoAmmend.txt", sep="")
  write.table(clustIn, file=inputFileName, sep=",", quote=FALSE, row.names=FALSE)
  
  CNVtotalList <- clustIn
  
  for(currCNV in 1:nrow(CNVtotalList)){
    if(CNVtotalList[currCNV, "type"] == "tetrasomy"){
      CNVtotalList[currCNV, "stats"] <- 2
      CNVtotalList[currCNV, "type"] <- 2
    }else if(CNVtotalList[currCNV, "type"] == "trisomy"){
      CNVtotalList[currCNV, "stats"] <- 1
      CNVtotalList[currCNV, "type"] <- 2
    }else if(CNVtotalList[currCNV, "type"] == "BLOH"){
      CNVtotalList[currCNV, "stats"] <- 0
      CNVtotalList[currCNV, "type"] <- 2
    }else{
      CNVtotalList[currCNV, "stats"] <- 0
      CNVtotalList[currCNV, "type"] <- 0
    }
  }
  CNVtotalList["start"] <- NULL
  CNVtotalList["end"] <- NULL
  CNVtotalList[10:11] <- CNVtotalList[6:7]
  CNVtotalList <- CNVtotalList[-c(1,6:7,12:14)]
  names(CNVtotalList) <- c("SetName", "CNV_ID", "Chromosome", "length", "a2", "b2", "alpha", "beta")
  
  #write timing file
  inputFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".matlab.input.txt", sep="")
  write.table(CNVtotalList, file=inputFileName, sep=",", quote=FALSE, row.names=FALSE)
  
  
}

