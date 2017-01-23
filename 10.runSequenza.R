#3.1.runSequenza.R
#
#
#
################### notes ###################
#
# call CNVs using sequenza
#
# 3.0.makeSequenzaScripts.R (3.1 and 3.2)
#           |
# run .sh scripts on apocrita (note: normal pileups have to be made first)
#           |
#           V
# 3.1.runSequenza.R
#           
#
#
################# libraries #################
library("sequenza")


################# subroutines #################



############### main program ################

#get sample list
#sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)

sampleNames <- unique(sampleList[[1]])

#dirs
sequenzaDir <- "22.sequenzaCalls/"
rawDir <- "rawFiles/"
seqFiles <- ".binned.seqz.gz"
GCplots <- "GCplots/"
anaDir <- "analysis/"

#set parameter search space
#cellTest <- seq(0.4,1,0.01)
#ploidyTest <- seq(1, 2, 0.1)

cellTest <- seq(0.6,1,0.01)
ploidyTest <- seq(1, 4, 0.1)

#make sequenza prep scripts (for each sample)
for(currSam in 1:nrow(sampleList)){
  print(paste("#### performing analysis for ", sampleList[currSam, 1], "-", sampleList[currSam, 2], " ####",sep=""))
  
  #current sample
  currName <- sampleList[currSam, 2]
  if(strsplit(currName, "_")[[1]][2] == "NORMAL" | strsplit(currName, "_")[[1]][2] == "N"){
    next
  }
  
  #current set
  currSet <- sampleList[currSam, 1]
  subSample <- sampleList[sampleList[[1]]==currSet, ]
  
  #read data from tumour seqz file name
  tumourSeq <- paste(sampleList[1,6], sequenzaDir, rawDir, currSet, "/", currName, seqFiles, sep="")
  if(!file.exists(tumourSeq)){
    next
  }
  #   seqz.data <- read.seqz(tumourSeq)
  #   
  #   #assess GC content
  #   gc.stats <- gc.norm(x = seqz.data$depth.ratio, gc = seqz.data$GC.percent)
  #   
  #   gc.vect <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
  #   seqz.data$adjusted.ratio <- seqz.data$depth.ratio / gc.vect[as.character(seqz.data$GC.percent)]
  #   
  #   #plot CG depth ratio corrections
  #   pdf(file=paste(sampleList[1,6], sequenzaDir, GCplots, currSet, "/", currName, ".GCplots.pdf", sep=""), height = 5, width = 10)
  #     par(mfrow = c(1,2), cex = 1, las = 1, bty ='l')
  #     matplot(gc.stats$gc.values, gc.stats$raw, type ='b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2), xlab ='GC content (%)', ylab ='Uncorrected depth ratio')
  #     legend('topright', legend = colnames(gc.stats$raw), pch = c(1, 19, 1))
  #     hist2(seqz.data$depth.ratio, seqz.data$adjusted.ratio, breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2), xlab ='Uncorrected depth ratio', ylab ='GC-adjusted depth ratio')
  #   dev.off()
  
  
  #analyze data
  system(command=paste("mkdir", paste(sampleList[1,6], sequenzaDir, anaDir, currSet, "/", currName, "/", sep="")))
  
  #remove headers included in table before analysis
  seqIn <- read.table(tumourSeq, header = TRUE, sep="\t", stringsAsFactors = FALSE)
  seqIn <- seqIn[seqIn[[1]]!="chromosome", ]
  tumourSeqOut <- paste(sampleList[1,6], sequenzaDir, rawDir, currSet, "/", currName, ".binned.seqz", sep="")
  write.table(seqIn, file = tumourSeqOut, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  system(command = paste("gzip -f", tumourSeqOut))

  test <- sequenza.extract(file = tumourSeq, chromosome.list = c(1:4,6:22, "X", "Y"), min.reads = 15)
  #seqzFile <- read.seqz(tumourSeq)
  
  #plot depth and BAF ratios (not needed, included below)
  #pdf(file=paste(sampleList[1,6], sequenzaDir, anaDir, currSet, "/", currName, "/", currName, ".BAF-depthRatio.pdf", sep=""), height = 5, width = 8)
  #  chromosome.view(mut.tab = test$mutations[[1]], baf.windows = test$BAF[[1]], ratio.windows = test$ratio[[1]], min.N.ratio = 1, segments = test$segments[[1]], main = test$chromosomes[1])
  #dev.off()
  
  #infer cellularity and ploidy
  CP.example <- sequenza.fit(test, cellularity = cellTest, ploidy = ploidyTest)
  
  #sequenza analysis
  sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = currName, out.dir=paste(sampleList[1,6], sequenzaDir, anaDir, "/", currSet, "/", currName, "/", sep="") )
  
}