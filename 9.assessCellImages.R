# takes sample list and cellCounts table as input
# performs linear regression on Ki67 and bCat expression

################### notes: #####################

#################### libraries #################

#################### subroutines #################

#################### main program #################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
setNames <- unique(sampleList[[1]])
setNames <- setNames[-c(12:24)]

#get cell count table
fileName <- "/Users/cross01/PhD/CRCproject/20.cellImageAnalysis/cellCounts.csv"
cellTab <- read.csv(file=fileName, stringsAsFactors = FALSE)

#calculate expression proportions
cellTab["Ki67prop"] <- (cellTab["Ki67"] + cellTab["Ki67_Negative"])
cellTab["Ki67prop"] <- cellTab["Ki67"] / cellTab["Ki67prop"] 

cellTab["B.catprop"] <- (cellTab["B.cat"] + cellTab["B.cat_Negative"])
cellTab["B.catprop"] <- cellTab["B.cat"] / cellTab["B.catprop"] 

cellTab[cellTab[["cancerRegion"]]!="CENTRE", "cancerRegion"] <- "OUTER"


#perfrom linear model of exonic burden and ploidy against ki67 and bcat proliferation

#ki67
ki67SNV <- lm(data=cellTab, cellTab[["exonicVariants"]] ~ cellTab[["Ki67prop"]])
ki67SNVSummary <- summary(ki67SNV)
ki67ploidy <- lm(data=cellTab, cellTab[["avPloidy"]] ~ cellTab[["Ki67prop"]])
ki67ploidySummary <- summary(ki67ploidy)

pdf(file="~/PhD/CRCproject/20.cellImageAnalysis/SNVvsKi67.pdf", width = 5, height = 5)
plot(x=cellTab[["Ki67prop"]], y=cellTab[["exonicVariants"]], ylab = "SNV exonic burden", xlab = "Ki67 %")
abline(ki67SNV, col="red")
text(x = 0.8, y= 82, labels = paste("pvalue:", ki67SNVSummary$coefficients[2, 4]) )
text(x = 0.8, y= 75, labels = paste("R2:", ki67SNVSummary$r.squared) )
dev.off()


pdf(file="~/PhD/CRCproject/20.cellImageAnalysis/ploidyvsKi67.pdf", width = 5, height = 5)
plot(x=cellTab[["Ki67prop"]], y=cellTab[["avPloidy"]], ylab = "ploidy", xlab = "Ki67 %")
abline(ki67ploidy, col="red")
text(x = 0.8, y= 2.5, labels = paste("pvalue:", ki67ploidySummary$coefficients[2, 4]) )
text(x = 0.8, y= 2.2, labels = paste("R2:", ki67ploidySummary$r.squared) )
dev.off()


#bcat
bcatSNV <- lm(data=cellTab, cellTab[["exonicVariants"]] ~ cellTab[["B.catprop"]])
bcatSNVSummary <- summary(bcatSNV)
bcatploidy <- lm(data=cellTab, cellTab[["avPloidy"]] ~ cellTab[["B.catprop"]])
bcatploidySummary <- summary(bcatploidy)


pdf(file="~/PhD/CRCproject/20.cellImageAnalysis/SNVvsbcat.pdf", width = 5, height = 5)
plot(x=cellTab[["B.catprop"]], y=cellTab[["exonicVariants"]], ylab = "SNV exonic burden", xlab = "Bcat %")
abline(bcatSNV, col="red")
text(x = 0.8, y= 82, labels = paste("pvalue:", bcatSNVSummary$coefficients[2, 4]) )
text(x = 0.8, y= 75, labels = paste("R2:", bcatSNVSummary$r.squared) )
dev.off()


pdf(file="~/PhD/CRCproject/20.cellImageAnalysis/ploidyvsbcat.pdf", width = 5, height = 5)
plot(x=cellTab[["B.catprop"]], y=cellTab[["avPloidy"]], ylab = "ploidy", xlab = "Bcat %")
abline(bcatploidy, col="red")
text(x = 0.8, y= 2.5, labels = paste("pvalue:", bcatploidySummary$coefficients[2, 4]) )
text(x = 0.8, y= 2.2, labels = paste("R2:", bcatploidySummary$r.squared) )
dev.off()




#stats table by sample comparisons
statsTab <- data.frame(matrix(0, nrow = length(setNames), ncol = 8))
names(statsTab) <- c("Ki67-cancerRegion", "Ki67-exonicVariants", "Ki67-MSIperc", "Ki67-avPloidy", "B.cat-cancerRegion", "B.cat-exonicVariants", "B.cat-MSIperc", "B.cat-avPloidy")
row.names(statsTab) <- setNames

for(currSam in 1:length(setNames)){
  
  #subset main table
  cellTemp <- cellTab[cellTab[["Set"]]==setNames[currSam], ]
  
  #test for adequate samples
  if(nrow(cellTemp) > 2){
    for(currTest in 1:ncol(statsTab)){
      #curr test variables
      varTemp <- names(statsTab)[currTest]
      var1 <- paste(strsplit(varTemp, split = "-")[[1]][1], "prop", sep="")
      var2 <- strsplit(varTemp, split = "-")[[1]][2]
        
      statsRes <- summary(lm(data=cellTemp, cellTemp[[var1]] ~ cellTemp[[var2]]))
      statsTab[currSam, currTest] <- paste(round(statsRes$coefficients[2,4], digits = 3), round(statsRes$adj.r.squared, digits = 3), sep=":")
    }
  }else{
    #add NAs to table and skip to next sample
    statsTab[currSam,] <- NA
    next()
  }
}

#save table
statsName <- "/Users/cross01/PhD/CRCproject/20.cellImageAnalysis/cellCounts.bySample.stats.csv"
write.csv(statsTab, file=statsName, quote = FALSE, row.names = TRUE)

cellTemp <- cellTab[cellTab[["cancerRegion"]]=="CENTRE", ]
cellTemp <- cellTemp[cellTemp[["Set"]]!="Set.04", ]

#stats table
statsTab <- data.frame(matrix(0, nrow = 2, ncol = 6))
names(statsTab) <- c("Ki67-exonicVariants", "Ki67-MSIperc", "Ki67-avPloidy", "B.cat-exonicVariants", "B.cat-MSIperc", "B.cat-avPloidy")
row.names(statsTab) <- c("CENTRE", "OUTER")

#stats comparing center biopsies only
statsTab["CENTRE", "Ki67-exonicVariants"] <- paste(round(summary(lm(data=cellTemp, Ki67prop ~ exonicVariants))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, Ki67prop ~ exonicVariants))$r.squared, digits = 3), sep=":")
statsTab["CENTRE", "Ki67-MSIperc"] <- paste(round(summary(lm(data=cellTemp, Ki67prop ~ MSIperc))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, Ki67prop ~ MSIperc))$r.squared, digits = 3), sep=":")
statsTab["CENTRE", "Ki67-avPloidy"] <- paste(round(summary(lm(data=cellTemp, Ki67prop ~ avPloidy))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, Ki67prop ~ avPloidy))$r.squared, digits = 3), sep=":")

statsTab["CENTRE", "B.cat-exonicVariants"] <- paste(round(summary(lm(data=cellTemp, B.catprop ~ exonicVariants))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, B.catprop ~ exonicVariants))$r.squared, digits = 3), sep=":")
statsTab["CENTRE", "B.cat-MSIperc"] <- paste(round(summary(lm(data=cellTemp, B.catprop ~ MSIperc))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, B.catprop ~ MSIperc))$r.squared, digits = 3), sep=":")
statsTab["CENTRE", "B.cat-avPloidy"] <- paste(round(summary(lm(data=cellTemp, B.catprop ~ avPloidy))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, B.catprop ~ avPloidy))$r.squared, digits = 3), sep=":")

cellTemp <- cellTab[cellTab[["cancerRegion"]]=="OUTER", ]
cellTemp <- cellTemp[cellTemp[["Set"]]!="Set.04", ]

statsTab["OUTER", "Ki67-exonicVariants"] <- paste(round(summary(lm(data=cellTemp, Ki67prop ~ exonicVariants))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, Ki67prop ~ exonicVariants))$r.squared, digits = 3), sep=":")
statsTab["OUTER", "Ki67-MSIperc"] <- paste(round(summary(lm(data=cellTemp, Ki67prop ~ MSIperc))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, Ki67prop ~ MSIperc))$r.squared, digits = 3), sep=":")
statsTab["OUTER", "Ki67-avPloidy"] <- paste(round(summary(lm(data=cellTemp, Ki67prop ~ avPloidy))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, Ki67prop ~ avPloidy))$r.squared, digits = 3), sep=":")

statsTab["OUTER", "B.cat-exonicVariants"] <- paste(round(summary(lm(data=cellTemp, B.catprop ~ exonicVariants))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, B.catprop ~ exonicVariants))$r.squared, digits = 3), sep=":")
statsTab["OUTER", "B.cat-MSIperc"] <- paste(round(summary(lm(data=cellTemp, B.catprop ~ MSIperc))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, B.catprop ~ MSIperc))$r.squared, digits = 3), sep=":")
statsTab["OUTER", "B.cat-avPloidy"] <- paste(round(summary(lm(data=cellTemp, B.catprop ~ avPloidy))[[4]][2,4], digits = 3), round(summary(lm(data=cellTemp, B.catprop ~ avPloidy))$r.squared, digits = 3), sep=":")

statsName <- "/Users/cross01/PhD/CRCproject/20.cellImageAnalysis/cellCounts.merged.stats.csv"
write.csv(statsTab, file=statsName, quote = FALSE, row.names = TRUE)



#analyse and plot region conparisons
dataIn <- read.table(file="~/PhD/CRCproject/17.imageAnalysis/summaryTable.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)
dataIn <- dataIn[c(1:73), ]
dataIn[dataIn[["cancerRegion"]]!="CENTRE", "cancerRegion"] <- "OUTER"
dataIn <- dataIn[dataIn[["setName"]]!="Set.04", ]


#outer vs centre SNV exon burden
tesRes <- wilcox.test(dataIn[dataIn[["cancerRegion"]]=="OUTER", "exonicVariants"], dataIn[dataIn[["cancerRegion"]]=="CENTRE", "exonicVariants"])

pdf(file="~/PhD/CRCproject/20.cellImageAnalysis/outerCentreSNVburden.pdf")
  boxplot(dataIn[dataIn[["cancerRegion"]]=="OUTER", "exonicVariants"], dataIn[dataIn[["cancerRegion"]]=="CENTRE", "exonicVariants"], names = c("OUTER", "CENTRE"))
  text(x = 1, y= 180, labels = tesRes$p.value)
dev.off()


#outer vs centre ploidy
tesRes <- wilcox.test(dataIn[dataIn[["cancerRegion"]]=="OUTER", "avPloidy"], dataIn[dataIn[["cancerRegion"]]=="CENTRE", "avPloidy"])

pdf(file="~/PhD/CRCproject/20.cellImageAnalysis/outerCentrePloidy.pdf")
  boxplot(dataIn[dataIn[["cancerRegion"]]=="OUTER", "avPloidy"], dataIn[dataIn[["cancerRegion"]]=="CENTRE", "avPloidy"], names = c("OUTER", "CENTRE"))
  text(x = 1, y= 3, labels = tesRes$p.value)
dev.off()


#SNV vs ploidy burdenmodel
modRes <- lm(dataIn[["exonicVariants"]] ~ dataIn[["avPloidy"]])
modSummary <- summary(modRes)
modSummary$coefficients[2, 4]

pdf(file="~/PhD/CRCproject/20.cellImageAnalysis/SNVvsPloidy.pdf")
  plot(y=dataIn[["exonicVariants"]], x=dataIn[["avPloidy"]], ylab = "SNV exonic burden", xlab = "ploidy")
  abline(modRes)
  text(x = 2.5, y= 180, labels = paste("pvalue:", modSummary$coefficients[2, 4]) )
  text(x = 2.5, y= 170, labels = paste("R2:", modSummary$r.squared) )
dev.off()


