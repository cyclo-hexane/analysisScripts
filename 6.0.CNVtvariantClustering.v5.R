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
    estBetaCenter <- estBetaCenter * celluEstSub 
    estAlphaCenter <- estAlphaCenter * celluEstSub 
    
    Afdata <- tempSNVtableSub[[5]]
    clusteredAF <- kmeans(Afdata, centers = c(estBetaCenter, estAlphaCenter))
    characterValues <- clusteredAF$cluster
    
    #stats uses  bootstrapping as per: http://stats.stackexchange.com/questions/11691/how-to-tell-if-data-is-clustered-enough-for-clustering-algorithms-to-produce-m
    bootstrapSub <- clusterboot(Afdata, B=100, bootmethod="boot", clustermethod=kmeansCBI, krange=2, centers = c(estBetaCenter, estAlphaCenter))
    bootstrapSub <- paste(round(bootstrapSub$bootmean, digits = 3), collapse = ":")
    
    #get cluster centre means and assign cluster charatcers to alpha or beta
    alphaValues <- Afdata[which(characterValues == 2)]
    betaValues <- Afdata[which(characterValues == 1)]
    
    alphaMean <- mean(alphaValues)
    betaMean <- mean(betaValues)
    
    distBeta <- abs(x = betaMean - estBetaCenter) 
    distAlpha <- abs(x = alphaMean - estAlphaCenter) 
    
    alpha <- length(alphaValues)
    beta <- length(betaValues)
    ciAlpha <- sd(alphaValues) * 1
    ciBeta <- sd(betaValues) * 1
    ciHighAlpha <- alphaMean + ciAlpha
    ciLowAlpha <- alphaMean - ciAlpha
    ciHighBeta <- betaMean + ciBeta
    ciLowBeta <- betaMean - ciBeta
    
    #get clustering partition point
    clusterPart <- max(betaValues)
    
    #correct alpha and beta using covSub
    binAlphaModel <- dbinom(x = ceiling(clusterPart * covSub):floor(covSub), size = round(covSub), prob = alphaMean)
    sdAlphaModel <- sd(x = binAlphaModel)
    alpha.portion <- sum(binAlphaModel)
    alpha.factor <- 1 / alpha.portion
    
    binBetaModel<- dbinom(ceiling(0.01 * covSub):floor(clusterPart * covSub), round(covSub), betaMean)
    sdBetaModel <- sd(x = binBetaModel)
    beta.portion <- sum(binBetaModel)
    beta.factor <- 1 / beta.portion
    
    alpha.as.beta.portion <- sum(dbinom(ceiling(0.01 * covSub):floor(clusterPart * covSub), round(covSub), alphaMean))
    beta.as.alpha.portion <- sum(dbinom(ceiling(clusterPart * covSub):floor(covSub), round(covSub), betaMean))
    beta.lost.portion <- sum(dbinom(0:floor(0.01 * covSub), round(covSub), betaMean))
    
    newAlpha <- round((alpha * alpha.factor) - (beta * beta.factor * beta.as.alpha.portion))
    newBeta <- round((beta * beta.factor) - (alpha * alpha.factor * alpha.as.beta.portion))
    
    #plot graph of clusters
    pdf(file=paste(graphOutDirSub, CNVstateSub, ".pdf", sep=""))
    par(xpd=FALSE)
    
      #plots
      hist(alphaValues, xlab="VAF", breaks=seq(0,1,0.01), xlim=c(0,1), ylim=c(0,15), col=rgb(1,0,0,0.5), main="clustering VAF distribution", freq = FALSE)
      denAlpha <- hist(alphaValues, breaks=seq(0,1,0.01), plot = FALSE)$density 
      hist(betaValues, xlab="VAF", breaks=seq(0,1,0.01), col=rgb(0,0,1,0.5), freq = FALSE, add=TRUE)
      denBeta <- hist(betaValues, breaks=seq(0,1,0.01), plot = FALSE)$density 
      
      #get parameters and plor overlayed gaussian model
      alphaDensity <- density(rnorm(n = 1000000, mean=alphaMean, sd=sdAlphaModel), kernel = "gaussian", from = -1, to = 1)
      maxAlphaYdensity <- max(alphaDensity$y)
      lines(alphaDensity, col="red3", pch=2, lwd=2)
      
      #get R2 distances and calculate R2 stats
      R2distances <- matrix(NA, nrow=length(denAlpha), ncol=4)
      R2distances[,1] <- denAlpha
      R2distances[,2] <- seq(0.01,1,0.01)
      for(currAdd in 1:nrow(R2distances)){
        R2distances[currAdd, 4] <- alphaDensity$x[findInterval(R2distances[currAdd, 2], vec = alphaDensity$x)]
        R2distances[currAdd, 3] <- alphaDensity$y[findInterval(R2distances[currAdd, 2], vec = alphaDensity$x)]
      }
      R2distances <- R2distances[R2distances[,1]>0, ]
      R2results <- summary(lm(R2distances[,3] ~ R2distances[,1]))
      
      betaDensity <- density(rnorm(n = 1000000, mean=betaMean, sd=sdBetaModel), kernel = "gaussian", from = -1, to = 1)
      maxBetaYdensity <- max(betaDensity$y)
      lines(betaDensity, col="blue4", pch=2, lwd=2)
      
      #get R2 distances and calculate R2 stats
      R2Betadistances <- matrix(NA, nrow=length(denBeta), ncol=4)
      R2Betadistances[,1] <- denBeta
      R2Betadistances[,2] <- seq(0.01,1,0.01)
      for(currAdd in 1:nrow(R2Betadistances)){
        R2Betadistances[currAdd, 4] <- betaDensity$x[findInterval(R2Betadistances[currAdd, 2], vec = betaDensity$x)]
        R2Betadistances[currAdd, 3] <- betaDensity$y[findInterval(R2Betadistances[currAdd, 2], vec = betaDensity$x)]
      }
      R2Betadistances <- R2Betadistances[R2Betadistances[,1]>0, ]
      R2Betaresults <- summary(lm(R2Betadistances[,3] ~ R2Betadistances[,1]))
      
      
      #save stats to table, catching linear regressions that cannot be performed
      if(nrow(R2results$coefficients) > 1){
        clusterStats <- paste(round(R2results$coefficients[2,4], digits = 4), ":", round(R2Betaresults$coefficients[2,4], digits = 4), sep="")
      }else{
        clusterStats <- NA
      }
      
      #plot cluster means
      lines(x=c(alphaMean, alphaMean), y=c(0, maxAlphaYdensity), lty=1, pch=1, lwd=4, col="red3")
      lines(x=c(betaMean, betaMean), y=c(0, maxBetaYdensity), lty=1, pch=1, lwd=4, col="blue4")
      
      #plot estimated cluster centres
      lines(x=c(estBetaCenter, estBetaCenter), y=c(0, 15), lty=2, pch=2, lwd=2, col="blue4")
      lines(x=c(estAlphaCenter, estAlphaCenter), y=c(0, 15), lty=2, pch=2, lwd=2, col="red3")
      
      #add CI for clustered dat
      lines(x = c(ciLowAlpha, ciHighAlpha), y=c(maxAlphaYdensity, maxAlphaYdensity), lty=1, pch=2, lwd=3, col="red3") 
      lines(x = c(ciLowBeta, ciHighBeta), y=c(maxBetaYdensity, maxBetaYdensity), lty=1, pch=2, lwd=3, col="blue4")
      
      #add cluster bootstrap values
      text(x=0.8, y=1, paste("R2 (a/b):", clusterStats), cex=0.6)
      
      #add cellularity values
      text(x=0.8, y=2, paste("cellularity:", celluEstSub), cex=0.6)
      
    dev.off()
    
    #check if means are within expected sd boundary
    if(distAlpha < ciAlpha & distBeta < ciBeta){
      #do nothing
    }else{
      clusterStats <- NA
    }
    
    #return parameters
    parameterList <- c(newAlpha, newBeta, bootstrapSub, clusterStats, CNVstateSub, alphaMean, betaMean, clusterPart, paste(graphOutDirSub, CNVstateSub, ".pdf", sep=""))
    
  }else{
    #cannot calculate small regions
    parameterList <- c(nrow(tempSNVtableSub), 0, NA, NA, CNVstateSub, NA, NA, NA, NA)
  }
  return(parameterList)
}

#### ammend CNVs with manually set VAF partitions
#tempSNVtableSub <- currVAFs; currentStateSub <- currStates; graphOutDirSub <- graphOutDir; partitionSub <- newPartition; covSub <- covData; celluEstSub<- cellularityEst
ammendParameters <- function(tempSNVtableSub, currentStateSub, graphOutDirSub, celluEstSub, partitionSub, covSub){
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
    estBetaCenter <- estBetaCenter * celluEstSub 
    estAlphaCenter <- estAlphaCenter * celluEstSub 
    
    Afdata <- tempSNVtableSub[[5]]
    
    #get cluster centre means and assign cluster charatcers to alpha or beta
    alphaValues <- Afdata[Afdata >= partitionSub]
    betaValues <- Afdata[Afdata < partitionSub]
    
    alphaMean <- mean(alphaValues)
    betaMean <- mean(betaValues)
    
    distBeta <- abs(x = betaMean - estBetaCenter) 
    distAlpha <- abs(x = alphaMean - estAlphaCenter) 
    
    alpha <- length(alphaValues)
    beta <- length(betaValues)
    ciAlpha <- sd(alphaValues) * 1
    ciBeta <- sd(betaValues) * 1
    ciHighAlpha <- alphaMean + ciAlpha
    ciLowAlpha <- alphaMean - ciAlpha
    ciHighBeta <- betaMean + ciBeta
    ciLowBeta <- betaMean - ciBeta
    
    #get clustering partition point
    clusterPart <- max(betaValues)
    
    #correct alpha and beta using covSub
    binAlphaModel <- dbinom(x = ceiling(clusterPart * covSub):floor(covSub), size = round(covSub), prob = alphaMean)
    sdAlphaModel <- sd(x = binAlphaModel)
    alpha.portion <- sum(binAlphaModel)
    alpha.factor <- 1 / alpha.portion
    
    binBetaModel<- dbinom(ceiling(0.01 * covSub):floor(clusterPart * covSub), round(covSub), betaMean)
    sdBetaModel <- sd(x = binBetaModel)
    beta.portion <- sum(binBetaModel)
    beta.factor <- 1 / beta.portion
    
    alpha.as.beta.portion <- sum(dbinom(ceiling(0.01 * covSub):floor(clusterPart * covSub), round(covSub), alphaMean))
    beta.as.alpha.portion <- sum(dbinom(ceiling(clusterPart * covSub):floor(covSub), round(covSub), betaMean))
    beta.lost.portion <- sum(dbinom(0:floor(0.01 * covSub), round(covSub), betaMean))
    
    newAlpha <- round((alpha * alpha.factor) - (beta * beta.factor * beta.as.alpha.portion))
    newBeta <- round((beta * beta.factor) - (alpha * alpha.factor * alpha.as.beta.portion))
    
    #plot graph of clusters
    pdf(file=paste(graphOutDirSub, CNVstateSub, ".pdf", sep=""))
    par(xpd=FALSE)
    
    #plots
    hist(alphaValues, xlab="VAF", breaks=seq(0,1,0.01), xlim=c(0,1), ylim=c(0,15), col=rgb(1,0,0,0.5), main="clustering VAF distribution", freq = FALSE)
    denAlpha <- hist(alphaValues, breaks=seq(0,1,0.01), plot = FALSE)$density 
    hist(betaValues, xlab="VAF", breaks=seq(0,1,0.01), col=rgb(0,0,1,0.5), freq = FALSE, add=TRUE)
    denBeta <- hist(betaValues, breaks=seq(0,1,0.01), plot = FALSE)$density 
    
    #get parameters and plor overlayed gaussian model
    alphaDensity <- density(rnorm(n = 1000000, mean=alphaMean, sd=sdAlphaModel), kernel = "gaussian", from = -1, to = 1)
    maxAlphaYdensity <- max(alphaDensity$y)
    lines(alphaDensity, col="red3", pch=2, lwd=2)
    
    #get R2 distances and calculate R2 stats
    R2distances <- matrix(NA, nrow=length(denAlpha), ncol=4)
    R2distances[,1] <- denAlpha
    R2distances[,2] <- seq(0.01,1,0.01)
    for(currAdd in 1:nrow(R2distances)){
      R2distances[currAdd, 4] <- alphaDensity$x[findInterval(R2distances[currAdd, 2], vec = alphaDensity$x)]
      R2distances[currAdd, 3] <- alphaDensity$y[findInterval(R2distances[currAdd, 2], vec = alphaDensity$x)]
    }
    R2distances <- R2distances[R2distances[,1]>0, ]
    R2results <- summary(lm(R2distances[,3] ~ R2distances[,1]))
    
    betaDensity <- density(rnorm(n = 1000000, mean=betaMean, sd=sdBetaModel), kernel = "gaussian", from = -1, to = 1)
    maxBetaYdensity <- max(betaDensity$y)
    lines(betaDensity, col="blue4", pch=2, lwd=2)
    
    #get R2 distances and calculate R2 stats
    R2Betadistances <- matrix(NA, nrow=length(denBeta), ncol=4)
    R2Betadistances[,1] <- denBeta
    R2Betadistances[,2] <- seq(0.01,1,0.01)
    for(currAdd in 1:nrow(R2Betadistances)){
      R2Betadistances[currAdd, 4] <- betaDensity$x[findInterval(R2Betadistances[currAdd, 2], vec = betaDensity$x)]
      R2Betadistances[currAdd, 3] <- betaDensity$y[findInterval(R2Betadistances[currAdd, 2], vec = betaDensity$x)]
    }
    R2Betadistances <- R2Betadistances[R2Betadistances[,1]>0, ]
    R2Betaresults <- summary(lm(R2Betadistances[,3] ~ R2Betadistances[,1]))
    
    
    #save stats to table, catching linear regressions that cannot be performed
    if(nrow(R2results$coefficients) > 1){
      clusterStats <- paste(round(R2results$coefficients[2,4], digits = 4), ":", round(R2Betaresults$coefficients[2,4], digits = 4), sep="")
    }else{
      clusterStats <- NA
    }
    
    #plot cluster means
    lines(x=c(alphaMean, alphaMean), y=c(0, maxAlphaYdensity), lty=1, pch=1, lwd=4, col="red3")
    lines(x=c(betaMean, betaMean), y=c(0, maxBetaYdensity), lty=1, pch=1, lwd=4, col="blue4")
    
    #plot estimated cluster centres
    lines(x=c(estBetaCenter, estBetaCenter), y=c(0, 15), lty=2, pch=2, lwd=2, col="blue4")
    lines(x=c(estAlphaCenter, estAlphaCenter), y=c(0, 15), lty=2, pch=2, lwd=2, col="red3")
    
    #add CI for clustered dat
    lines(x = c(ciLowAlpha, ciHighAlpha), y=c(maxAlphaYdensity, maxAlphaYdensity), lty=1, pch=2, lwd=3, col="red3") 
    lines(x = c(ciLowBeta, ciHighBeta), y=c(maxBetaYdensity, maxBetaYdensity), lty=1, pch=2, lwd=3, col="blue4")
    
    #add cluster bootstrap values
    text(x=0.8, y=1, paste("R2 (a/b):", clusterStats), cex=0.6)
    
    #add cellularity values
    text(x=0.8, y=2, paste("cellularity:", celluEstSub), cex=0.6)
    
    dev.off()
    
    #check if means are within expected sd boundary
    if(distAlpha < ciAlpha & distBeta < ciBeta){
      #do nothing
    }else{
      clusterStats <- NA
    }
    
    #return parameters
    parameterList <- c(newAlpha, newBeta, NA, clusterStats, CNVstateSub, alphaMean, betaMean, clusterPart, paste(graphOutDirSub, CNVstateSub, ".pdf", sep=""))
    
  }else{
    #cannot calculate small regions
    parameterList <- c(nrow(tempSNVtableSub), 0, NA, NA, CNVstateSub, NA, NA, NA, NA)
  }
  return(parameterList)
}


################## libraries ##################

library(fpc)
library(mixtools)

################## main program ##################



#global variables
# 1. sampleList file
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.Ha.01.filt.csv", header=FALSE, stringsAsFactors=FALSE)

# 2. CNV segmentation file
CNVPrepended <- ".penalty0.95.baf.gt.txt"

# 3. SNV copy state file
#variantCallFile <- ".snv.annoVar.variant_function.txt"
variantCallFile <- ".snv.annoVar.variant_function.0.01.txt"

# 4. purity estimates
purityFile <- ".penalty0.95.purity-cna-baf.txt"
workingDir <- "7.CNVcalls.archive/baf.updated/"
#SNVworkingDir <- "1.platypusCalls/somaticTotal/"
SNVworkingDir <- "1.platypusCalls/somaticTotal.0.01/"
analysisDir <- "6.CNVtiming/inputData.final/"
purityDir <- "13.purityEstimates/"
setNames <- unique(sampleList[[1]])
setNames <- setNames[c(2,3,5,9,10,11,16)]

# segmentations to exclude
excDir <- "6.CNVtiming/exclusionRegions/"
excludeSegs <- c(7)

#vcf list for ammendments
vcfList <- as.list(NA)

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
  samplesTimed <- sampleNames
  
  #get cellularity estimates
  celluEst <- read.table(file=paste(currentSamples[1,6], purityDir, currentSamples[1, 1], purityFile, sep=""), sep=" ")
  
  remSamples <- celluEst[celluEst[[4]]<=0.5, 3]
  if(length(remSamples) > 0){
    for(currExc in 1:length(remSamples)){
      samplesTimed <- samplesTimed[samplesTimed != remSamples[currExc]]
    }
  }
  
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
  
  #add line if needed
  #vcfData <- as.data.frame(append(vcfData, values = NA, after = 0))
  
  #remove normal sample
  vcfData <- vcfData[-c((8 + normalSample) , (8 + (noSamples + 1) + normalSample))]
  vcfData[1] <- c(1:nrow(vcfData))
  
  #calculate allele frequencies and change X/Y to 23/24
  for(y in 1:noSamples){
    vcfData[(ncol(vcfData)+1)] <- as.vector(vcfData[(8+(noSamples)+y)]) / as.vector(vcfData[8+y])
  }
  names(vcfData) <- c("line", "region", "variant", "chrom", "pos", "pos2", "ref", "alt", paste(sampleNames, ".NR", sep=""), paste(sampleNames, ".NV", sep=""), paste(sampleNames, ".AF", sep=""))
  chromData <- vcfData[["chrom"]]
  #chromData[chromData=="X"] <- 23
  #chromData[chromData=="Y"] <- 24
  vcfData["chrom"] <- chromData
  
  #remove variants not in CNV regions that are hypermutated regions, set to blank if there are no regions to remove
  if(excludeSegs == currSet){
    excFile <- paste(currentSamples[1,6], excDir, currentSamples[1,1], ".txt", sep="")
    excTab <- read.table(file=excFile, sep="\t", stringsAsFactors = FALSE, header=FALSE)
  }else{
    excTab <- data.frame(0, nrow=0, ncol=3)
  }
  
  vcfDataUnfilt <- vcfData
  
  #remove snv if not truncal or in exclusion table above
  vcfCutoff <- 0
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
    vcfTabCutoff <- table(vcfData[currSNV, paste(sampleNames, ".AF", sep="")] > vcfCutoff)
    if(!("TRUE" %in% names(vcfTabCutoff))){
      remRow[remCounter] <- currSNV
      remCounter <- remCounter + 1
      next()
    }else if(as.numeric(vcfTabCutoff["TRUE"]) != noSamples){
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
  
  vcfList[[currSet]] <- vcfData
  
  
  
  ####### C. loop through segmentations and attempt to cluster ####### 
  
  #setup results table
  CNVtotalList <- data.frame(matrix(NA, nrow=(length(samplesTimed) * nrow(filteredSegCNVs)), ncol=18))
  names(CNVtotalList) <- c("SetName", "SampleID", "CNV.ID", "Chromosome", "CNVlength", "start", "end", "alpha", "beta", "bootstraps", "stats", "type", "alphaMean", "betaMean", "clusterPartition", "phyloLoc", "cellularity", "filter")
  CNVtotalList[1] <- currentSamples[1,1]
  CNVtotalList[2] <- samplesTimed

  CNVcounter <- 1
  colSeqMaj <- paste(samplesTimed, "_Major", sep="")
  colSeqMin <- paste(samplesTimed, "_Minor", sep="")
  for(currCNV in 1:nrow(filteredSegCNVs)){
    currChrom <- filteredSegCNVs[currCNV, "chr"]
    currPosStart <- filteredSegCNVs[currCNV, "first.locus"]
    currPosEnd <-  filteredSegCNVs[currCNV, "last.locus"]
    CNVlength <- currPosEnd - currPosStart
    
    #for each sample with CNV present get SNVs, calculate alpha and beta and add to table
    for(currSam in 1:length(samplesTimed)){
      currMajCol <- colSeqMaj[currSam]
      currMinCol <- colSeqMin[currSam]
      currentState <- filteredSegCNVs[currCNV, c(currMinCol, currMajCol)]
      if(filteredSegCNVs[currCNV, "chr"] == 23){
        intervalID <- paste("X:", (filteredSegCNVs[currCNV, "first.locus"]+1), "-", filteredSegCNVs[currCNV, "last.locus"], sep="")
      }else{
        intervalID <- paste(filteredSegCNVs[currCNV, "chr"], ":", (filteredSegCNVs[currCNV, "first.locus"]+1), "-", filteredSegCNVs[currCNV, "last.locus"], sep="")
      }
      
      currID <- paste(samplesTimed[currSam], "-chr", currChrom, "-", currentState[1], ".", currentState[2], "-", CNVlength, sep="")
      
      print(paste("##### processing ID:", currID, "#####") )
      
      #current cellularity value
      cellularityEst <- celluEst[celluEst[[3]]==samplesTimed[currSam], 4]
      
      # 1. get SNV for each region
      tempSNVtable <- getSNVs(currPosStart, currPosEnd, currChrom, samplesTimed[currSam], vcfData)
      
      #get coverage information from biopsy specific file
      covFileName <- paste('~/PhD/CRCproject/6.CNVtiming/gatk_coverage/', setNames[currSet], '_', samplesTimed[currSam], '.sample_interval_summary', sep = '')
      covData <- read.table(covFileName, header = T, fill = T)
      covData <- covData[covData[["Target"]]==intervalID, "average_coverage"]
      
      # 2. calculate alpha and beta
      graphOutDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currID, ".", sep="")
      parameters <- calcParameters(tempSNVtable, currentState, graphOutDir, cellularityEst, covData)
      
      #if parameters not possible NA is placed instead of calculated values for alpha and beta
      CNVtotalList[CNVcounter,] <- c(currentSamples[1,1], samplesTimed[currSam], currID, currChrom, CNVlength, currPosStart, currPosEnd, parameters[1:8], filteredSegCNVs[currCNV, "phyloLoc"], cellularityEst, "NA")
      
      #check for filter criteria and mark accordingly, also move graph to correct folder
      filterString <- c()
      filtCounter <- 1
      if(nrow(tempSNVtable) < 50){
        filterString <- "lowVariants"
      }else if(is.na(parameters[4])){
        #mean of alpha cluster out of confidence 
        filterString <- "clustering"
        newDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/exclude/", sep="") 
        system(command = paste("mv", parameters[9], newDir))
      }else{
        #pvaluesClust <- as.numeric(strsplit(parameters[4], ":")[[1]])
        #if(pvaluesClust[1] > 0.05 | pvaluesClust[2] > 0.05){
        #  filterString[filtCounter] <- "clustering"
        #  filtCounter <- filtCounter + 1
        #}
        if(parameters[5] == "polysomy" | parameters[5] == "loss"){
          filterString[filtCounter] <- "copyState"
          filtCounter <- filtCounter + 1
        }
        if(filtCounter != 1){
          newDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/exclude/", sep="") 
          system(command = paste("mv", parameters[9], newDir))
        }else{
          filterString[1] <- "PASS"
          newDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/include/", sep="")
          system(command = paste("mv", parameters[9], newDir))
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
    CNVtotalList <- CNVtotalList[c(2:5,9:10,6:7)]
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
  




####### D. check segmentation plots manually and move segmentations not passed to include folder ####### 

#loop though and reassess include graphs
for(currSet in 2:length(setNames)){
  currentSamples <- subset(sampleList, sampleList[1]==setNames[currSet])
  
  celluEst <- read.table(file=paste(currentSamples[1,6], purityDir, currentSamples[1, 1], purityFile, sep=""), sep=" ")
  
  print(paste("re-assessing set", currentSamples[1,1]))
  
  inputClustFile <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".input.parameters.total.txt", sep="")
  clustIn <- read.table(file=inputClustFile, sep=",", header = TRUE, stringsAsFactors = FALSE)
  
  #get include event list from graphs in include directory
  includeDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/include/", sep="")
  graphsIncluded <- list.files(path = includeDir, all.files = FALSE)
  
  #parse names
  tempFileHold <- matrix(unlist(strsplit(graphsIncluded, "\\.")), ncol = 4, nrow = length(graphsIncluded), byrow = TRUE)
  CNVsIncluded <- paste(tempFileHold[,1], tempFileHold[,2], sep=".")
  
  for(currCNV in 1:length(CNVsIncluded)){
    clustIn[clustIn[["CNV.ID"]]==CNVsIncluded[currCNV], "filter"] <- "INC"
  }
  
  #loop through ammended lists with manually stipulated partitions
  ammendFile <- paste(currentSamples[1,6], "6.CNVtiming/ammendLists/", currentSamples[1,1], ".ammend.txt", sep="")
  if(file.exists(ammendFile)){
    #get ammend list
    ammendData <- read.table(file=ammendFile, sep="\t", header = FALSE, stringsAsFactors = FALSE)
    
    #get total CNV timing list (to get alpha and beta values)
    totalFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".input.parameters.total.txt", sep="")
    totalData <- read.table(file=totalFileName, sep=",", header = TRUE, stringsAsFactors = FALSE)
    
    #get truncal vcf table
    currVCF <- vcfList[[currSet]]
    
    #loop through each CNV and ammend alpha and beta by the manually selected partition
    for(currAmm in 2:nrow(ammendData)){
      #get inputs for timing model
      currID <- ammendData[currAmm, 1]
      currSetName <- strsplit(currID, split = "-")[[1]][1]
      
      newPartition <- ammendData[currAmm, 2]
      currSampleName <- paste(currSetName, ".AF", sep="")
      currState <- strsplit(currID, split = "-")[[1]][3]
      currStates <- strsplit(currState, "\\.")[[1]]
      graphOutDir <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currID, ".", sep="")
      
      #get coverage datas
      covFileName <- paste('~/PhD/CRCproject/6.CNVtiming/gatk_coverage/', setNames[currSet], '_', currSetName, '.sample_interval_summary', sep = '')
      covData <- read.table(covFileName, header = T, fill = T)
      intervalTemp <- totalData[totalData[["CNV.ID"]]==currID, c("Chromosome", "start", "end", "CNVlength")]
      intervalID <- paste(intervalTemp[1,1], ":", (intervalTemp[1,2]+1), "-", intervalTemp[1,3], sep="")
      covData <- covData[covData[["Target"]]==intervalID, "average_coverage"]
      
      #get cellularity estimate
      cellularityEst <- celluEst[celluEst[[3]]==currSetName, 4]
      
      #get vcfs for this region
      currVAFs <- getSNVs(intervalTemp[1,2], intervalTemp[1,3], intervalTemp[1,1], currSetName, currVCF)
      
      parameters <- ammendParameters(currVAFs, currStates, graphOutDir, cellularityEst, newPartition, covData)
      
      clustIn <- rbind(clustIn, c(currentSamples[1,1], currSetName, currID, intervalTemp[1,1], intervalTemp[1,4], intervalTemp[1,2], intervalTemp[1,3], parameters[1:8], "T", cellularityEst, "INC"))
    }
  }
  
  
  
  #filter for variants with PASS
  clustIn <- clustIn[clustIn[["filter"]]=="INC",]
  clustIn["filter"] <- NULL
  inputFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".input.parameters.ammend.txt", sep="")
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
  CNVtotalList <- CNVtotalList[c(2:5,9:10,6:7)]
  names(CNVtotalList) <- c("SetName", "CNV_ID", "Chromosome", "length", "a2", "b2", "alpha", "beta")
  
  #write timing file
  inputFileName <- paste(currentSamples[1,6], analysisDir, currentSamples[1,1], "/", currentSamples[1,1], ".matlab.input.ammend.txt", sep="")
  write.table(CNVtotalList, file=inputFileName, sep=",", quote=FALSE, row.names=FALSE)
  
  
}

