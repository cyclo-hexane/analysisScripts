# script takes sampleList and cloneHD segmentations as input
# outputs table of phylogenetically and ancestrally correct CNV changes
# and CNV loads for rate and barcharts


#################    notes   ################# 
# the essense of this method is to use the simplest CNV transitions to explain the observed segmentation states
#
#     main program
#           |
#     for each chromosome               ----------- if simple score without model
#           |                           |                                                                 -------------- 2.getCNVcombinations
#     1.call modelCNVs <--------------->                                                                  |
#           |                           |---------- if complex and < 6 segs, score with parsimony model ---------------- 3.getParsimonyScore
#     merge all tables                  |                                                                 |
#                                       |                                                                 -------------- 4.assessModels
#                                       |---------- if complex and > 6 segs, score by classification
#                                                                                                   
#                                                                                                           


#################  libraries  #################


#################  functions  ################# 

# for each given CNV event get all possible CNV orders
#CNVlistSub<-CNVlist; setNamesTemp <- setNamesSub
getCNVcombinations <- function(CNVlistSub, setNamesTemp){
  CNVlistSub["ID"] <- paste("ID", CNVlistSub[["ID"]], sep="")
  
  setTemp <- strsplit(setNamesTemp[1], split="_")[[1]][1]
  print(paste("getting combinations for sample ", setTemp, "chromosome", CNVlistSub[1, "chr"]))
  noSamplesTemp <- length(setNamesTemp)
  noEvents <- nrow(CNVlistSub)
  noCNVtypes <- unique(CNVlistSub[c("Minor", "Major")])
  
  #add frequency to table
  noCNVtypes[3] <- NA
  names(noCNVtypes)[3] <- "freq"
  for(currType in 1:nrow(noCNVtypes)){
    noCNVtypes[currType, "freq"] <- nrow(CNVlistSub[CNVlistSub["Minor"]==noCNVtypes[currType, "Minor"] & CNVlistSub["Major"]==noCNVtypes[currType, "Major"], ])
  }
  
  #add combinable events for each CNV (treated separately)
  l <- rep(list(1:noEvents), noEvents)
  allComb <- expand.grid(l)
  names(allComb) <- paste("ID", c(1:nrow(CNVlistSub)), sep="")
  
  #now remove those not valid for this dataset
  keepRow <- rep(NA, nrow(allComb))
  for(testRow in 1:nrow(allComb)){
    assessIDs <- allComb[testRow,]
    assessRow <- table(c(1:noEvents) %in% allComb[testRow, ])
    if(as.numeric(assessRow["TRUE"])==noEvents){
      keepRow[testRow] <- TRUE 
    }else{
      #test to see whether this is a viable combinatory event
      testNumber <- unique(as.numeric(assessIDs))
      testNumber <- c(1:length(testNumber))
      testNumberOrders <- table(as.numeric(assessIDs))
      testNumberOrders <- testNumberOrders[order(as.numeric(names(testNumberOrders)))]
      if(identical(as.numeric(testNumber), as.numeric(names(testNumberOrders)))){
        #if number orders are correct
        validVect <- matrix(NA, nrow=1, ncol=noEvents)
        colnames(validVect) <- CNVlistSub[["ID"]]
        for(currValid in 1:length(testNumber)){
          currSegNames <- colnames(assessIDs)[assessIDs==as.numeric(currValid)]
          currRes <- CNVlistSub[CNVlistSub[["ID"]] %in% currSegNames, ]
          if(length(unique(currRes[["Minor"]]))==1 & length(unique(currRes[["Major"]]))==1){
            #if all regions are the same state mark as valid for this step
            validVect[assessIDs==currValid] <- 1
          }else{
            validVect[assessIDs==currValid] <- 0
          }
        }
        if(sum(validVect)==noEvents){
          keepRow[testRow] <- TRUE
        }else{
          keepRow[testRow] <- FALSE
        }
      }else{
        keepRow[testRow] <- FALSE
        next
      }  
    }
  }
  allComb <- allComb[keepRow, ] 
  nCombRows <-nrow(allComb)
  
  #identify 'normal' CNV segmentations
  normalSegs <- c()
  norCounter <- 1
  for(currAden in 1:nrow(CNVlistSub)){
    currMajor <- CNVlistSub[currAden, "Major"]
    currMinor <- CNVlistSub[currAden, "Minor"]
    if(currMajor==1 & currMinor==1){
      normalSegs[norCounter] <- currAden
      norCounter <- norCounter + 1
    }
  }
  
  #remove models where normals are gained first
  if(!is.null(normalSegs)){
    for(remNor in 1:length(normalSegs)){
      allComb <- allComb[!allComb[[normalSegs[remNor]]]==1, ]
    }
  }
  
  #make output table
  names(allComb) <- paste("ID", c(1:noEvents), sep="")  
  row.names(allComb) <- paste("Model_", c(1:nrow(allComb)), sep="")  
  return(allComb)
}



#for a given order assess parsimony score
#CNVcombTabSub <- CNVcombTab; CNVlistSub <- CNVlist; setNamesTemp <- setNamesSub; finalCheckStates <- assessSeg
getParsimonyScore <- function(CNVcombTabSub, CNVlistSub, setNamesTemp, finalCheckStates){
  CNVlistSub["ID"] <- paste("ID", c(1:ncol(CNVcombTabSub)), sep="")
  
  #print progress
  setTemp <- strsplit(setNamesTemp[1], split="_")[[1]][1]
  print(paste("getting parsimony scores for sample ", setTemp, "chromosome", CNVlistSub[1, "chr"]))
  noSamplesTemp <- length(setNamesTemp)
  
  #add score column
  CNVcombTabSub[ncol(CNVcombTabSub)+1] <- NA
  names(CNVcombTabSub)[ncol(CNVcombTabSub)] <- "score"
  
  CNVcombTabSub[ncol(CNVcombTabSub)+1] <- NA
  names(CNVcombTabSub)[ncol(CNVcombTabSub)] <- "stateScore"
  
  CNVlistSubCopy <- CNVlistSub
  
  CNVtabList <- as.list(NA)
  
  #for each model get parsimony score and check validity against finalCheckStates
  for(currMod in 1:nrow(CNVcombTabSub)){
    #setup segmentation region states table
    tempStates  <- finalCheckStates
    tempStates[5:ncol(tempStates)] <- 1
    parsScore <- 0
    stateScore <- 0
    CNVlistSub <- CNVlistSubCopy
    CNVlistSub[setNames] <- 0
    delRow <- c()
    delCounMain <- 1
    
    #iterate through model and calculated score 
    for(currChg in 1:(ncol(CNVcombTabSub)-1) ){
      #check if change exists
      if(!currChg %in% CNVcombTabSub[currMod,]){
        next
      }
      newStates <- tempStates
      
      #get change states
      currSegChanges <- names(CNVcombTabSub)[which(CNVcombTabSub[currMod, ]==currChg)]
      
      #get and enact current change
      changeSamples <- CNVlistSubCopy[CNVlistSub[["ID"]] %in% currSegChanges, ]
      minEnactStates <- changeSamples[1, "Minor"]
      majEnactStates <- changeSamples[1, "Major"]
      startEnactRegion <- as.numeric(min(changeSamples["first.locus"]))
      endEnactRegion <- as.numeric(max(changeSamples["last.locus"])) 
      
      #if >1 segment is changed ensure encapsulated regions are changed as well  
      if(nrow(changeSamples) > 1){
        #get samples effected
        tempNameList <- as.list(NA)
        maxNoSamples <- 0
        representCheck <- matrix(0, nrow=1, ncol=length(setNamesTemp))
        colnames(representCheck) <- setNamesTemp
        
        #now check which samples are effected by these regions
        for(currAdd in 1:nrow(CNVlistSubCopy)){
          firstLociTest <- as.numeric(CNVlistSubCopy[currAdd, "first.locus"])
          lastLociTest <- as.numeric(CNVlistSubCopy[currAdd, "last.locus"])
          testStates <- newStates[newStates[["first.locus"]] >= firstLociTest & newStates[["last.locus"]] <= lastLociTest, ]
          testfinalStates <- finalCheckStates[finalCheckStates[["first.locus"]] >= firstLociTest & finalCheckStates[["last.locus"]] <= lastLociTest, ]
          
          #if there is another segmentation in the list that is within this region
          if(firstLociTest >= startEnactRegion & lastLociTest <= endEnactRegion){
            #and this region is not normal (diploid)
            if(CNVlistSubCopy[currAdd, "Major"] != 1 | CNVlistSubCopy[currAdd, "Minor"] == 0){
              addSam <- setNames[CNVlistSubCopy[currAdd, setNames] == 1]
              #check samples are not at final state, do not add to addSam if so
              samRem <- c()
              remCounter <- 1
              for(currStateAss in 1:length(addSam)){
                testSamStates <- testStates[ c(paste(addSam[currStateAss], "_Minor", sep=""), paste(addSam[currStateAss], "_Major", sep="")) ]
                testSamStatesFinal <- testfinalStates[ c(paste(addSam[currStateAss], "_Minor", sep=""), paste(addSam[currStateAss], "_Major", sep="")) ]
                if(FALSE %in% (testSamStates == testSamStatesFinal)){
                  #do not remove from change list
                }else{
                  samRem[remCounter] <- currStateAss 
                  remCounter <- remCounter + 1
                }
              }
              #is there any samples to remove from the add list?
              if(!is.null(samRem)){
                addSam <- addSam[-samRem]
              }
              if(!length(addSam)==0){
                representCheck[1, addSam] <- 1
              }             
            }
          }
        }
        if(sum(representCheck)!=0){
          #column names
          firstCols <- paste(setNames[representCheck[1,]==1], "_Minor", sep="")
          lastCols <- paste(setNames[representCheck[1,]==1], "_Major", sep="") 
          
          #assess number of samples to be changed
          maxNoSamples <- as.numeric(table(representCheck==1)["TRUE"])
          
          #assign phylogenetic location
          if(maxNoSamples==noSamplesTemp){
            phyloDes <- "T"
          }else if(maxNoSamples==1){
            phyloDes <- "L"
          }else{
            phyloDes <- "B"
          }
          
          #mark new states in relivent regions
          firstLociEnactStates <- min(changeSamples["first.locus"])
          lastLociEnactStates <- max(changeSamples["last.locus"])
          
          #check each region in turn and change as needed
          for(currEnact in 1:nrow(newStates)){
            currSegStart <- newStates[currEnact, "first.locus"]
            currSegFin <- newStates[currEnact, "last.locus"]
            if(currSegStart >= firstLociEnactStates & currSegFin <= lastLociEnactStates){
              newStates[currEnact, firstCols] <- minEnactStates
              newStates[currEnact, lastCols] <- majEnactStates
            }
          }
        }else{
          #there were no changes to make
        }
      }else{
        #if only one region to change
        representCheck <- matrix(0, nrow=1, ncol=length(setNamesTemp))
        colnames(representCheck) <- setNamesTemp
        
        #mark samples to be changed
        representCheck[1, setNamesTemp[changeSamples[setNamesTemp]==1]] <- 1
        
        #add samples with non-normal states in same segmentation region, but only if not already at final state
        for(currAdd in 1:nrow(CNVlistSubCopy)){
          firstLociTest <- as.numeric(CNVlistSubCopy[currAdd, "first.locus"])
          lastLociTest <- as.numeric(CNVlistSubCopy[currAdd, "last.locus"])
          testStates <- newStates[newStates[["first.locus"]] >= firstLociTest & newStates[["last.locus"]] <= lastLociTest, ]
          testfinalStates <- finalCheckStates[finalCheckStates[["first.locus"]] >= firstLociTest & finalCheckStates[["last.locus"]] <= lastLociTest, ]
          
          #if there is another segmentation in the list that is within this region
          if(firstLociTest >= as.numeric(changeSamples["first.locus"]) & lastLociTest <= as.numeric(changeSamples["last.locus"])){
            #and this region is not normal (diploid)
            if(CNVlistSubCopy[currAdd, "Major"] != 1 | CNVlistSubCopy[currAdd, "Minor"] == 0){
              addSam <- setNames[CNVlistSubCopy[currAdd, setNames] == 1]
              #check samples are not at final state, do not add to addSam if so
              samRem <- c()
              remCounter <- 1
              for(currStateAss in 1:length(addSam)){
                testSamStates <- testStates[ c(paste(addSam[currStateAss], "_Minor", sep=""), paste(addSam[currStateAss], "_Major", sep="")) ]
                testSamStatesFinal <- testfinalStates[ c(paste(addSam[currStateAss], "_Minor", sep=""), paste(addSam[currStateAss], "_Major", sep="")) ]
                if(FALSE %in% (testSamStates==testSamStatesFinal)){
                  #do nothing
                }else{
                  samRem[remCounter] <- currStateAss 
                  remCounter <- remCounter + 1
                }
              }
              #is there any samples to remove from the add list?
              if(!is.null(samRem)){
                addSam <- addSam[-samRem]
              }
              if(!length(addSam)==0){
                representCheck[1, addSam] <- 1
              }             
            }
          }
        }
        
        #column names
        firstCols <- paste(setNames[representCheck[1,]==1], "_Minor", sep="")
        lastCols <- paste(setNames[representCheck[1,]==1], "_Major", sep="") 
        
        #get phylogenetic locations
        maxNoSamples <- as.numeric(table(representCheck==1)["TRUE"])
        if(maxNoSamples==noSamplesTemp){
          phyloDes <- "T"
        }else if(maxNoSamples==1){
          phyloDes <- "L"
        }else{
          phyloDes <- "B"
        }
        
        firstLociEnactStates <- changeSamples[1, "first.locus"]
        lastLociEnactStates <- changeSamples[1, "last.locus"]
        
        #check each region in turn and change as needed
        for(currEnact in 1:nrow(newStates)){
          currSegStart <- newStates[currEnact, "first.locus"]
          currSegFin <- newStates[currEnact, "last.locus"]
          if(currSegStart >= firstLociEnactStates & currSegFin <= lastLociEnactStates){
            newStates[currEnact, firstCols] <- minEnactStates
            newStates[currEnact, lastCols] <- majEnactStates
          }
        }
      }
      
      #compare new states to old and add to parsimony score 
      changeRecord <- sum(abs(newStates - tempStates))
      stateScore <- stateScore + changeRecord
      if(changeRecord >= 1){
        parsScore <- parsScore + 1
      }
      
      tempStates <- newStates
      
      #if no changes are needed, seg is ignored in this model
      if(changeRecord==0){
        phyloDes <- "removed"
      }
      
      #mark new phylogenetic locationa and effected samples
      CNVlistSub[CNVlistSub[["ID"]] %in% currSegChanges, "phyloLoc"] <- phyloDes
      representCheck <- representCheck[1, representCheck > 0]
      CNVlistSub[CNVlistSub[["ID"]] %in% currSegChanges, names(representCheck)] <- 1
      
      #merge multiple transitions to one row
      if(nrow(changeSamples) > 1){
        maxChg <- max(CNVlistSub[CNVlistSub[["ID"]] %in% currSegChanges, "last.locus"])
        CNVlistSub[CNVlistSub[["ID"]] %in% currSegChanges[1], "last.locus"] <- maxChg
        currSegChanges <- currSegChanges[currSegChanges!=currSegChanges[1]]
        CNVlistSub[CNVlistSub[["ID"]] %in% currSegChanges, "phyloLoc"] <- "removed"
      } 
    }
    
    #check validity of final state, mark model for removal if not
    if(FALSE %in% (finalCheckStates == tempStates)){
      #mark for deletion
      CNVcombTabSub[currMod, "score"] <- 0
      CNVcombTabSub[currMod, "stateScore"] <- 0
    }else{
      #add score to row
      CNVcombTabSub[currMod, "score"] <- parsScore
      CNVcombTabSub[currMod, "stateScore"] <- stateScore 
    }
        
    #save CNVlistSub to list
    CNVtabList[[currMod]] <- CNVlistSub
    names(CNVtabList)[currMod] <- row.names(CNVcombTabSub)[currMod]
  }
  
  returnList <- as.list(NA)
  returnList[[1]] <- CNVcombTabSub
  returnList[[2]] <- CNVtabList
  return(returnList)
}



#for parsimony model list, get best model and correct CNVlist table
#parseModelSub <- modelProp;CNVlistSub <- modelProposals; setNamesTemp <- setNamesSub; dirOutSub <- holdingSub
assessModels <- function(parseModelSub, CNVlistSub, setNamesTemp, dirOutSub){
  #criteria:
  # 1. assume most phylogenetically 'correct' order to be favorable i.e seg order should be Truncal -> Branch -> Leaf
  # 2. parsimony score (PS) and state change score (SCS) are minimized (SCS = individual changes are made per model-step)
  # 3. correct models for 'removed' segmentations
  # 4. prefer first state changes to be small
  # 5. check that best models are phylogenetically the same, if so choice is irrelevant
  # 6. after all else take fist from list, print error if more than once solution
  
  noSegSub <- ncol(parseModelSub)-2
  setTemp <- strsplit(setNamesTemp[1], split="_")[[1]][1]
  noSamplesSub <- length(setNamesTemp)
  chromosomeName <- CNVlistSub[[1]][1, "chr"]
  
  # 1. filter models by phylogenetic criteria
  modelNames <- paste(rownames(parseModelSub), "_phyloLoc", sep="")
  modeList <- rep(NA, nrow(parseModelSub))
  names(modeList) <- rownames(parseModelSub)
  phyloTreeOrder <- c("T", "B", "L")
  for(currPhy in 1:nrow(parseModelSub)){
    phyloTestOrder <- CNVlistSub[[row.names(parseModelSub)[currPhy]]][order(parseModelSub[currPhy, c(1:noSegSub)]), "phyloLoc"]
    phyloOrder <- phyloTestOrder[phyloTestOrder!="removed"]
    phyloOrder <- match(phyloOrder, phyloTreeOrder)
    if(all(phyloOrder == cummax(phyloOrder))){
      modeList[currPhy] <- 1
    }else{
      modeList[currPhy] <- 0
    } 
  }
  if(1 %in% modeList){
    parseModelSub <- parseModelSub[modeList==1, ]
  }else{
    print(paste("cannot filter models by phylogenetic order for", setTemp, "chr", CNVlistSub[1,"chr"],": filter ignored"))
    #write.table(minParseModel, file=paste(dirOutSub, ".fittingModels.txt", sep=""))  
  }
  
  
  # 2. filter models by state change (SCS) and parsimony score
  minPS <- min(parseModelSub["score"])
  minParseModel <- parseModelSub[parseModelSub[["score"]]==minPS, ]
  if(nrow(minParseModel)!=0){
    minSCS <- min(minParseModel["stateScore"])
    minSCSmodel <- minParseModel[minParseModel[["stateScore"]]==minSCS, ]
    if(nrow(minSCSmodel)!=0){
      parseModelSub <- minSCSmodel
    }else{
      print(paste("cannot filter models by PS for", setTemp, "chr", CNVlistSub[1,"chr"],": filter ignored"))
      parseModelSub <- minParseModel
    }
  }else{
    print(paste("cannot filter models by SCS for", setTemp, "chr", CNVlistSub[1,"chr"],": filter ignored"))
  }
  #remove score columns
  parseModelSub["score"] <- NULL
  parseModelSub["stateScore"] <- NULL
  
  
  
  # 3. correct models for removed events
  for(currCorr in 1:nrow(parseModelSub)){
    assessRow <- CNVlistSub[[row.names(parseModelSub)[currCorr]]]
    if("removed" %in% assessRow[["phyloLoc"]]){
      #model contains a 'removed' segmentation
      removedID <- assessRow[which(assessRow["phyloLoc"]=="removed"), "ID"]
      parseModelSub[currCorr, removedID] <- NA
      nonRem <- which(!names(parseModelSub) %in% removedID)
      nonRem <- parseModelSub[currCorr, nonRem]
      for(currAss in 1:length(unique(as.numeric(nonRem)))){
        if(TRUE %in% (nonRem == currAss)){
          next
        }else{
          testNo <- nonRem[nonRem > currAss]
          testNo <- min(testNo)
          nonRem[which(nonRem %in% testNo)] <- currAss
          parseModelSub[currCorr, which(parseModelSub[currCorr,] %in% testNo)] <- currAss
        }
      }
    }
  }
  parseModelSub <- unique(parseModelSub)
  
  
  # 4. prefer first state changes to be small i.e 1,3 is preferred to 2,3 as a first state change
  stateChangeList <- c()
  SCcounter <- 1
  for(currSC in 1:nrow(parseModelSub)){
    tempChange <- names(parseModelSub)[parseModelSub[currSC,]==1]
    tempChange <- tempChange[!is.na(tempChange)]
    currModName <- row.names(parseModelSub)[currSC]
    stateTemp <- CNVlistSub[[currModName]][which(CNVlistSub[[currModName]][["ID"]]%in%tempChange), c("Minor", "Major")]
    stateChangeList[SCcounter] <- abs(sum(unique(stateTemp)) - 2)
    SCcounter <- SCcounter + 1
  }
  if(length(min(stateChangeList))!=1){
    parseModelSub <- parseModelSub[which(stateChangeList==min(stateChangeList)), ]
  }
  
  
  # 5. check if best models are all phylogenetically the same
  if(nrow(parseModelSub) != 1){
    phyloTable <- parseModelSub
    for(currApp in 1:nrow(phyloTable)){
      phyloTable[currApp, ] <- CNVlistSub[[row.names(phyloTable)[currApp]]][["phyloLoc"]]
    }
    if(nrow(unique(phyloTable)) == 1){
      parseModelSub <- parseModelSub[1, ]
    }
  }
  
    
  # 6. take fist as best model (if models are phylogenetically the same this doesn't matter)
  if(nrow(parseModelSub) != 1){
    print(paste("######## cannot find best model for", setTemp, "chr", CNVlistSub[[1]][1,"chr"],"given filters ########"))
    write.table(parseModelSub, file=paste(dirOutSub, "/bestFittingModels.txt", sep=""), sep="\t", quote = FALSE, col.names = TRUE)
    bestModel <- parseModelSub[1,]
  }else{

    bestModel <- parseModelSub[1,]
  }
  
  #now correct CNVlist with model
  bestModelName <- rownames(bestModel)
  CNVlistSub <- CNVlistSub[[bestModelName]]
  
  #populate with modelled CNVs
  bestModel[is.na(bestModel)] <- "X"
  
  #assess whether to merge segmentations, remove those not needed
  if("removed" %in% CNVlistSub[["phyloLoc"]]){
    remRow <- which("removed" == CNVlistSub[["phyloLoc"]])
    delRow <- c()
    delCounter <- 1
    for(currMergeRows in 1:length(remRow)){
      upperRow <- NULL
      lowerRow <- NULL
      
      #if segmentations exist above and below to seg begin removed, test for merging
      if(remRow[currMergeRows] > 1){
        upperRow <- CNVlistSub[remRow[currMergeRows]-1, ]
      }
      if(remRow[currMergeRows] < nrow(CNVlistSub)){
        lowerRow <- CNVlistSub[remRow[currMergeRows]+1, ]
      }
      if(!is.null(upperRow) & !is.null(lowerRow)){
        if(identical(as.numeric(upperRow[c("Minor", "Major", setNames)]), as.numeric(lowerRow[c("Minor", "Major", setNames)]))){
          #merge upper and lower segmentations if states match and ether interval is the same or state is the same
          if(identical(as.numeric(CNVlistSub[remRow[currMergeRows], c("first.locus", "last.locus")]), as.numeric(upperRow[c("first.locus", "last.locus")])) |  identical(as.numeric(CNVlistSub[remRow[currMergeRows], c("first.locus", "last.locus")]), as.numeric(lowerRow[c("first.locus", "last.locus")])) | identical(as.numeric(CNVlistSub[remRow[currMergeRows], c("Minor", "Major")]), as.numeric(upperRow[c("Minor", "Major")])) ){
            CNVlistSub[remRow[currMergeRows]-1, "last.locus"] <- CNVlistSub[remRow[currMergeRows]+1, "last.locus"]
            delRow[delCounter] <- remRow[currMergeRows]+1
            delCounter <- delCounter + 1
          }
        }
      }
    }
    CNVlistFinal <- CNVlistSub[-c(remRow, delRow), ]
    CNVlistFinal[9+noSamplesSub] <- bestModelName
  }else{
    #assess whether segmentations can be merged
    delRow <- c()
    delCounter <- 1
    for(currMergeRows in 2:nrow(CNVlistSub)){
      if(identical(as.numeric(CNVlistSub[currMergeRows, c("Minor", "Major", setNames)]), as.numeric(CNVlistSub[(currMergeRows-1), c("Minor", "Major", setNames)])) ){
        CNVlistSub[(currMergeRows-1), "last.locus"] <- CNVlistSub[currMergeRows, "last.locus"]
        delRow[delCounter] <- currMergeRows
        delCounter <- delCounter + 1
      }
    }
    if(!is.null(delRow)){
      CNVlistSub <- CNVlistSub[-delRow, ]  
    }
    CNVlistFinal <- CNVlistSub
  }  
  #correct names columns
  CNVlistFinal[9+noSamplesSub] <- bestModelName
  names(CNVlistFinal)[9+noSamplesSub] <- "model"
  names(CNVlistFinal)[8+noSamplesSub] <- "phyloLoc"
  
  #write best model
  bestModelName <- paste(dirOutSub, "/bestModel.chr", chromosomeName,".scoreTab.txt", sep="")
  write.table(bestModel, file=bestModelName, sep="\t", row.names=TRUE, col.names = TRUE, quote = FALSE)
  
  #return CNVlist table
  return(CNVlistFinal)
}



#main CNV model functions
#assessSeg <- segSub; setNamesSub <- setNames; noSamplesSub <- noSamples; holdingSub <- tabDir
modelCNVs <- function(assessSeg, setNamesSub, noSamplesSub, holdingSub){
  #### notes ####
  # function takes a list of segmentations and produces an ancestrally and parsimoniously derived CNV table
  # 1. coverts to 'sample-wise' table of events with binary states
  # 2. if simple, i.e one CNV interval and CNV gain, assign phylogenetic state as is
  # 3. if not, use parsimony method to determine most likely order
  
  noSamplesSub <- length(setNamesSub)
  noSegmentations <- nrow(assessSeg)
  sampleSet <- strsplit(setNamesSub, "_")
  sampleSet <- sampleSet[[1]][1]
  chromName <- assessSeg[1, "chr"]
  
  # 1. separate segmentations to get recovered events for each sample
  CNVlist <- data.frame(matrix(NA, nrow=10, ncol=(8+noSamplesSub)))
  names(CNVlist) <- c("ID", "chr", "first.locus", "nloci", "last.locus", "Minor", "Major", setNamesSub, "phyloLoc")
  IDcounter <- 1
  for(currRow in 1:nrow(assessSeg)){
    assessMajRow <- assessSeg[currRow, paste(setNamesSub, "_Major", sep="")]
    assessMinRow <- assessSeg[currRow, paste(setNamesSub, "_Minor", sep="")]
    combState <- paste(assessMinRow, ":", assessMajRow, sep="")
    for(addStates in unique(combState)){
      tempStates <- strsplit(addStates, ":")
      sampleStates <- setNamesSub[combState==addStates]
      CNVlist[IDcounter, "ID"] <- IDcounter
      CNVlist[IDcounter, "chr"] <- assessSeg[currRow, 1]
      CNVlist[IDcounter, "Minor"] <- as.numeric(tempStates[[1]][1])
      CNVlist[IDcounter, "Major"] <- as.numeric(tempStates[[1]][2])
      CNVlist[IDcounter, "first.locus"] <- assessSeg[currRow, "first.locus"]
      CNVlist[IDcounter, "last.locus"] <- assessSeg[currRow, "last.locus"]
      CNVlist[IDcounter, "nloci"] <- (CNVlist[IDcounter, "last.locus"] - CNVlist[IDcounter, "first.locus"]) / 1000000
      CNVlist[IDcounter, sampleStates] <- 1
      if(length(sampleStates)!=noSamplesSub){
        zeroStates <- setNamesSub[!setNamesSub %in% sampleStates]
        CNVlist[IDcounter, zeroStates] <- 0
      }
      IDcounter <- IDcounter + 1
    }
  }
  CNVlist <- CNVlist[complete.cases(CNVlist[-ncol(CNVlist)]), ]
  
  ######### try to merge segmentations, if possible #########
  delRow <- c()
  delCounter <- 1
  CNVstates <- unique(CNVlist[,c("Minor", "Major")])
  
  #use Bp differences between intervals as test for matching
  regionIntervals <- unique(assessSeg[,c("first.locus", "last.locus")])
  regionDiffs <- data.frame(matrix(NA, ncol=2, nrow=(nrow(regionIntervals)-1)))
  for(currDiff in 1:nrow(regionDiffs)){
    regionDiffs[currDiff, 1] <- regionIntervals[(currDiff+1), 1] - regionIntervals[currDiff, 1] 
    regionDiffs[currDiff, 2] <- regionIntervals[(currDiff+1), 2] - regionIntervals[currDiff, 2]
  }
  #now test the segmentation differences for each state
  for(currDel in 1:nrow(CNVstates)){
    #subset by segs with this state
    currSegAssess <- CNVlist[CNVlist[["Minor"]]==CNVstates[currDel, "Minor"] & CNVlist[["Major"]]==CNVstates[currDel, "Major"], ]
    
    #do segmentations effect the same sample?
    if(nrow(unique(currSegAssess[setNames]))!=1){
      orderFlag <- 1
      next
    }
    
    #assess whether segmentations are sequential
    assIntervals <- currSegAssess[c("first.locus", "last.locus")]
    if(nrow(assIntervals)==1){
      next
    }
    orderFlag <- 0
    for(currAssIn in 1:(nrow(assIntervals)-1)){
      #which segmentation is it from regionIntervals 
      seg1 <- which(assIntervals[currAssIn, 1] %in% regionIntervals[[1]])
      seg2 <- which(assIntervals[currAssIn, 2] %in% regionIntervals[[2]])
      if(seg1 != seg2){
        #if intervals are not the same, this region cannot be merged
        orderFlag <- 1
        next
      }else{
        currDiffInt <- assIntervals[(currAssIn+1), ] - assIntervals[currAssIn, ]
        if(!identical(as.numeric(currDiffInt), as.numeric(regionDiffs[currAssIn, ]))){
          #if difference intervals do not match regionDiffs
          orderFlag <- 1
        }
      }
    }
    
    #merge region or not?
    if(orderFlag == 1){
      #cannot merge region
      next
    }else{
      #merge segmentations for current state
      CNVlist[CNVlist["ID"]==currSegAssess[1, "ID"], "last.locus"] <- CNVlist[CNVlist["ID"]==currSegAssess[nrow(currSegAssess), "ID"], "last.locus"]
      
      #mark those to be deleted
      for(delAs in 2:nrow(currSegAssess)){
        delRow[delCounter] <- currSegAssess[delAs, "ID"]
        delCounter <- delCounter + 1
      }
    }
  }
  #delete rows by ID
  if(!is.null(delRow)){
    CNVlist <- CNVlist[-delRow, ]
    
    #rename IDs
    CNVlist["ID"] <- c(1:nrow(CNVlist))
  }
  
  ######### 2. if the CNV is 'simple' #########
  
  #remove normal segs and test for 'simplicity'
  remList <- c()
  remCounter <- 1
  for(currRem in 1:nrow(CNVlist)){
    remAssess <- CNVlist[currRem, c("Minor" , "Major", setNamesSub)]
    if(as.numeric(remAssess["Minor"])==1 & as.numeric(remAssess["Major"])==1){
      remList[remCounter] <- currRem
      remCounter <- remCounter + 1
    }
  }
  if(length(remList) > 0){
    testCNV <- CNVlist[-remList, ]
  }else{
    testCNV <- CNVlist
  }
  #assign phylogenetic regions using frequency alone
  if(length(unique(testCNV[["first.locus"]]))==1 & length(unique(testCNV[["last.locus"]]))==1 & length(unique(testCNV[["Minor"]]))==1 & length(unique(testCNV[["Major"]]))==1){
    samplesPresent <- setNamesSub[testCNV[setNamesSub]==1]
    if(length(samplesPresent)==noSamplesSub){
      testCNV[1, "phyloLoc"] <- "T"
    }else if(length(samplesPresent)==1){
      testCNV[1, "phyloLoc"] <- "L"
    }else{
      testCNV[1, "phyloLoc"] <- "B"
    }
    testCNV[1, samplesPresent] <- 1
    CNVlist <- testCNV
  }else{
    CNVlist["phyloLoc"] <- "overlapping"
  }
  
  
  ######### 3. if events are non-simple, use parsimony to determine most likely ancestral states #########
  if("overlapping" %in% CNVlist[["phyloLoc"]] & nrow(CNVlist) < 8 & nrow(CNVlist) > 1){
    #tabulate all CNV combinations
    CNVcombTab <- getCNVcombinations(CNVlist, setNamesSub)
    
    
    #get parsimony scores
    scoreTab <- getParsimonyScore(CNVcombTab, CNVlist, setNamesSub, assessSeg)
    
    #remove non-valid models
    modelProp <- scoreTab[[1]]
    modelProp <- modelProp[modelProp[["stateScore"]]!=0, ]
    
    #make directory to hold data
    system(command=paste("mkdir ", holdingSub, sep=""))
    
    #output model list and phylogenetic locations
    phyloTab <- paste(holdingSub, "/allModels.chr", chromName,".scoreTab.txt", sep="")
    write.table(modelProp, file=phyloTab, sep="\t", quote = FALSE, row.names=TRUE, col.names = TRUE)
    modNamesList <- row.names(modelProp)
    
    #output model proposals
    modelProposals <- scoreTab[[2]]
    for(currModName in modNamesList){
      currModelProp <- modelProposals[[currModName]]
      modTabRes <- paste(holdingSub, "/", currModName,".modResTab.txt", sep="")
      write.table(currModelProp, file=modTabRes, sep="\t", quote = FALSE, row.names=FALSE, col.names = TRUE)
    }
    
    #determine best model and correct CNVlist table with result
    finalParseModel <- assessModels(modelProp, modelProposals, setNamesSub, holdingSub)
    
  }else if( ("overlapping" %in% CNVlist[["phyloLoc"]] & nrow(CNVlist) > 7) | ("overlapping" %in% CNVlist[["phyloLoc"]] & nrow(CNVlist)==1) ){
    #model using basic classifier
    remListSub <- c()
    remListCounter <- 1
    for(currClass in 1:nrow(CNVlist)){
      assessData <- CNVlist[currClass, c("Minor", "Major")]
      assessRow <- table(CNVlist[currClass, setNamesSub]==1)
      if(as.numeric(assessData[1])==1 & as.numeric(assessData[2])==1){
        #mark for removal
        remListSub[remListCounter] <- currClass
        remListCounter <- remListCounter + 1
      }else{
        #mark phylogenetic location
        if(as.numeric(assessRow["TRUE"])==noSamplesSub){
          CNVlist[currClass, "phyloLoc"] <- "T"
        }else if(as.numeric(assessRow["TRUE"])==1){
          CNVlist[currClass, "phyloLoc"] <- "L"
        }else{
          CNVlist[currClass, "phyloLoc"] <- "B"
        }
      }
    }
    #remove normal rows
    if(length(remListSub)!=0){
      finalParseModel <- CNVlist[-remListSub,]
    }else{
      finalParseModel <- CNVlist
    }
    
    #add "cannot model" column
    if(nrow(finalParseModel)==1){
      finalParseModel[ncol(finalParseModel)+1] <- "noModelNeeded"
    }else{
      finalParseModel[ncol(finalParseModel)+1] <- "cannotModel"
    }
    names(finalParseModel)[ncol(finalParseModel)] <- "model"
  }else{
    #do nothing no modelling required
    finalParseModel <- CNVlist
  }
  
  #correct columns (if not modelled)
  modName <- names(finalParseModel)[9+noSamples]
  if(is.na(modName)){
    finalParseModel[9+noSamplesSub] <- NA
    names(finalParseModel)[9+noSamplesSub] <- "model"
    finalParseModel[9+noSamplesSub] <- "noModelNeeded"
  }
  
  #return phylogenetic table
  finalParseModel <- finalParseModel[order(finalParseModel["first.locus"]), ]
  return(finalParseModel)
}  



#plotting function for inferred phylogenetic states
plotCNVstates <- function(mergedCNVTableSub, plotFileNameSub, chromoLengthSub){
  #remove 'y' chromosome from analysis
  chromoLengthSub <- chromoLengthSub[chromoLengthSub[1] != 24,]
  
  #make graph
  pdf(file=plotFileNameSub, 16, 2.5)
  par(mar=c(3,3,3,3))
  
  plot(-1, -1, xlim=c(0, 2839313410), ylim=c(0,5), xlab="chromosome", ylab="CNV state", xaxt="n", cex.axis=1.2, main="")
  for(currPhylo in c("T", "B", "L")){
    #get colour
    if(currPhylo=="T"){
      colPlot <- "steelblue"
    }else if(currPhylo=="B"){
      colPlot <- "goldenrod"
    }else{
      colPlot <- "salmon"
    }
    
    #plot each chromosome segment
    for(chr in c(1:23)){
      chrStart <- chromoLengthSub[chromoLengthSub[1]==chr, 3]
      
      #subset to get chromosome segments
      currChrom <- mergedCNVTableSub[mergedCNVTableSub$chr==chr & mergedCNVTableSub$phyloLoc==currPhylo, ]
      currChrom <- currChrom[complete.cases(currChrom),]
      if(nrow(currChrom)==0){
        next
      }
      for(segments in 1:nrow(currChrom)){
        #if overlapped set line to dotted 
        if(table(currChrom[currChrom[2]==currChrom[segments,2], 2]) > 1 & table(currChrom[currChrom[4]==currChrom[segments,4], 4]) > 1){
          overlap <- 2
          lineWeight <- 2
        }else{
          overlap <- 1
          lineWeight <- 4
        }
        
        ypoints <- as.numeric(currChrom[segments, "Minor"]) + as.numeric(currChrom[segments, "Major"]) + offSet
        segments(x0=chrStart+currChrom[segments, "first.locus"], x1=chrStart+currChrom[segments, "last.locus"], y0=ypoints, y1=ypoints, lty=overlap, lwd=lineWeight, col=colPlot)
      }
    }
  }  
  axis(1, at=chromoLengthSub[[4]],labels=c(1:23),cex.axis=0.8)
  abline(h=2, lty=2, col="grey", lwd=1)
  for(currLine in 2:nrow(chromoLengthSub)){
    abline(v=c(chromoLengthSub[currLine, 3]), lty=1, col="black", lwd=1)
  }
  dev.off()
}





################# main program ################# 


#sampleList <- read.csv(file="~/PhD/CRCproject/archive/masterSampleList.filtered.csv", header=FALSE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/PhD/CRCproject/8.CNVphylogenetics/trainingSampleList.csv", header=FALSE, stringsAsFactors=FALSE)

sampleNames <- unique(sampleList[[1]])

holdingDir <- "7.CNVcalls.final/baf.updated/"
fileApp <- ".penalty0.95.baf.gt.txt"
plotDir <- "8.CNVphylogenetics/"

#chromosome lengths file
chromoLength <- read.csv(file="~/PhD/ReferenceGenome/chromosomeStartEnd.csv", header=TRUE, colClasses=c(rep("numeric")))

for(currSamp in 1:length(sampleNames)){
  print(paste("######### starting analysis for sample set", sampleNames[currSamp], "##########"))
  
  #subset sampleList
  currSet <- sampleNames[currSamp]
  subSample <- subset(sampleList, sampleList[1]==currSet)
  noSamples <- subSample[1,8]-1
  normalIndex <- subSample[1,7]+1 
  setNames <- subSample[[2]]
  
  #remove normal from setNames
  normalName <- setNames[normalIndex]
  norRem <- which(normalName == setNames)
  setNames <- setNames[-norRem]
  subSample <- subSample[-norRem,]
  
  #get segmentation files
  segName <- paste(subSample[1,6], holdingDir, currSet, fileApp, sep="")
  segFile <- read.table(file=segName, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  segFile[segFile == "X"] <- 1

  #get segmentation length
  segFile["nloci"] <- segFile[["last.locus"]] - segFile[["first.locus"]]
  segFile["nloci"] <- segFile["nloci"] / 1000000
  
  #remove small segmentations by smoothing (currently set to 50 SNPs per seg)
  #focalMuts <- segFile[segFile["nloci"]<=9 & segFile["nloci"] > 9, ]
  segFile <- segFile[segFile["nloci"]>10, ]
  
  #make chromosome phylo table (to summerize phylogenetic loads)
  chrAbberations <- unique(segFile[[1]])
  chrAbTable <- data.frame(matrix(NA, ncol=13, nrow=length(chrAbberations)))
  names(chrAbTable) <- c("chr", "trunk", "branch", "leaf", "nTrunk", "nBranch", "nLeaf", "nLossTrunk", "nLossBranch", "nLossLeaf", "sizeLossTrunk", "sizeLossBranch", "sizeLossLeaf")
  chrAbTable[1] <- chrAbberations
  
  #make merge table for CNVs
  mergedCNVTable <- data.frame(matrix(NA, ncol=9+noSamples, nrow=0))
  names(mergedCNVTable) <- c("ID", "chr", "first.locus", "nloci", "last.locus", "Minor", "Major", setNames, "phyloLoc", "model")
  
  #assess load on trunk, branch, leaf for each chromosome (saving run info to log file)
  logFileName <- paste(subSample[1,6], plotDir, currSet, ".CNVmodel.log.txt", sep="")
  file.create(file1 = logFileName, showWarnings = FALSE)
  
  #assess load on trunk, branch, leaf for each chromosome
  for(currChrom in 1:nrow(chrAbTable)){
    print(paste("###### starting CNV model for", currSet, "chromosome", currChrom, "######"))
    
    currChr <- chrAbTable[currChrom, 1]
    segSub <- subset(segFile, segFile[1]==currChr)
    currLength <- chromoLength[chromoLength[[1]]==currChr, 2]
    
    #check if only normal segmentations
    normalFlag <- 0
    for(currCheck in 1:nrow(segSub)){
      assessCheck <- segSub[currCheck, ]
      rowTotals <- rep(NA, length(setNames))
      for(currSam in 1:length(setNames)){
        minCheck <- paste(setNames[currSam], "_Minor", sep="")
        majCheck <- paste(setNames[currSam], "_Major", sep="")
        if(assessCheck[1, minCheck]==1 & assessCheck[1, majCheck]==1){
          rowTotals[currSam] <- 0
        }else{
          rowTotals[currSam] <- 1
        }
      }
      if(sum(rowTotals) >= 1){
        normalFlag <- "T"
      }else if(sum(rowTotals) == 0 & normalFlag != "T" ){
        normalFlag <- "N"
      }
    }
    
    if(normalFlag == "N"){
      next
    }
    
    #model CNV changes (logging result)
    sink(file=logFileName, append=TRUE)
    tabDir <- paste(subSample[1,6], plotDir, currSet, "/", currSet, ".chr", currChrom, sep="")
    mergedTable <- modelCNVs(segSub, setNames, noSamples, tabDir)
    sink()
    
    #remove unwanted rows
    mergedTable <- mergedTable[mergedTable[["phyloLoc"]]!="removed", ]
    
    #merge CNV table
    if(nrow(mergedTable)!=0){
      mergedCNVTable <- rbind(mergedCNVTable, mergedTable)
      
      #total up phylogenetic loads for plotting
      trunkCNVs <- mergedTable[mergedTable["phyloLoc"]=="T",]
      branchCNVs <- mergedTable[mergedTable["phyloLoc"]=="B",]
      leafCNVs <- mergedTable[mergedTable["phyloLoc"]=="L",]
      noLossTrunk <- trunkCNVs[trunkCNVs[["Minor"]]==0, ]
      noLossBranch <- branchCNVs[branchCNVs[["Minor"]]==0, ]
      noLossLeaf <- leafCNVs[leafCNVs[["Minor"]]==0, ]
      chrAbTable[currChrom, c(2,5,8,11)] <- c(sum(trunkCNVs[["nloci"]]), nrow(trunkCNVs), nrow(noLossTrunk),  sum(noLossTrunk[["nloci"]]) )
      chrAbTable[currChrom, c(3,6,9,12)] <- c(sum(branchCNVs[["nloci"]]), nrow(branchCNVs), nrow(noLossBranch), sum(noLossBranch[["nloci"]]) )
      chrAbTable[currChrom, c(4,7,10,13)] <- c(sum(leafCNVs[["nloci"]]), nrow(leafCNVs), nrow(noLossLeaf), sum(noLossLeaf[["nloci"]]) )
    }
  }
  
  mergedCNVTable <- mergedCNVTable[order(mergedCNVTable["chr"]),]
  
  #write table
  mergedOut <- paste(subSample[1,6], plotDir, subSample[1,1], ".phyloCNVs.v2.csv" ,sep="")
  write.table(mergedCNVTable, file=mergedOut, sep=",", quote = FALSE, col.names = TRUE, row.names=FALSE)
  
  #write loads table
  loadOut <- paste(subSample[1,6], plotDir, subSample[1,1], ".CNVsLoads.csv" ,sep="")
  chrAbTable <- chrAbTable[complete.cases(chrAbTable), ]
  write.table(chrAbTable, file=loadOut, sep=",", quote = FALSE, col.names = TRUE, row.names=FALSE)
  
  
  #order CNV table
  #mergedCNVTable <- mergedCNVTable[order(mergedCNVTable["phyloLoc"]),]
  
  #plot CNV evolution to check above table
  #plotFileName <- paste(subSample[1,6], plotDir, subSample[1,1], ".phyloSegs.pdf" ,sep="")
  #plotCNVstates(mergedCNVTable, plotFileName, chromoLength)
}


