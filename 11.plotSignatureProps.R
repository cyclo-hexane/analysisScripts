# 11.plotSignatureProps.R


################# notes #################
# 1. read in EMu results files and plot mutational signature proportions
#
# 2. create stats for trunk-branch-leaf comparisons
#
#
#
################# main program #################

sampleList <- read.csv(file="~/PhD/CRCproject/masterSampleList.allSamples.filt.csv", header=FALSE, stringsAsFactors=FALSE)
sampleNames <- unique(sampleList[1])

holdingDir <- "5.mutationalSignatures/"
holdingDirVCF <- "1.platypusCalls/somaticTotal.0.01/"



################## 1. plot mutational signature proportions ################

#read in propotions file
propEMuRoot <- read.table(file=paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-Z.root.txt", sep=""), stringsAsFactors = FALSE, sep="\t")
propEMubranch <- read.table(file=paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-Z.branches.txt", sep=""), stringsAsFactors = FALSE, sep="\t")
propEMuleaf <- read.table(file=paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-Z.leaf.txt", sep=""), stringsAsFactors = FALSE, sep="\t")


#plot root signatures
pdfName <- paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-plot.root.pdf", sep="")
pdf(file=pdfName, onefile=TRUE, width=7, height=5)
  par(mar=c(7,5,5,5))
  plotMat <- t(as.matrix(propEMuRoot[4:7]))
  barplot(plotMat, col=c("olivedrab", "salmon", "royalblue", "goldenrod"), names.arg = propEMuRoot[[1]], las=2)
dev.off()


pdfName <- paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-plot.branches.pdf", sep="")
pdf(file=pdfName, onefile=TRUE, width=7, height=5)
par(mar=c(7,5,5,5))
plotMat <- t(as.matrix(propEMubranch[4:7]))
barplot(plotMat, col=c("olivedrab", "salmon", "royalblue", "goldenrod"), names.arg = propEMubranch[[1]], las=2)
dev.off()


pdfName <- paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-plot.leafs.pdf", sep="")
pdf(file=pdfName, onefile=TRUE, width=40, height=5)
par(mar=c(7,5,5,5))
plotMat <- t(as.matrix(propEMuleaf[4:7]))
barplot(plotMat, col=c("olivedrab", "salmon", "royalblue", "goldenrod"), names.arg = propEMuleaf[[1]], las=2)
dev.off()



#boxplots and stats of trunk-branch
rootComp <- propEMuRoot[-c(17:20, 24), c(4:7)]
names(rootComp) <- c("A", "B", "C", "D")
row.names(rootComp) <- propEMubranch[[1]]

branchComp <- propEMubranch[, c(4:7)]
names(branchComp) <- c("Ab", "Bb", "Cb", "Db")
row.names(branchComp) <- propEMubranch[[1]]

plotTab <- cbind(branchComp, rootComp)
plotTab <- plotTab[order(names(plotTab))]

#stats for carcinomas
carcinomasTab <- plotTab[c(1:3,5:11),]

sigAcomparison <- wilcox.test(carcinomasTab[["A"]], carcinomasTab[["Ab"]])$p.value
sigBcomparison <- wilcox.test(carcinomasTab[["B"]], carcinomasTab[["Bb"]])$p.value
sigCcomparison <- wilcox.test(carcinomasTab[["C"]], carcinomasTab[["Cb"]])$p.value
sigDcomparison <- wilcox.test(carcinomasTab[["D"]], carcinomasTab[["Db"]])$p.value

pdfName <- paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-boxplots.carcinomas.pdf", sep="")
pdf(file=pdfName, onefile=TRUE, width=5, height=5)
par(mar=c(7,5,5,5), xpd=TRUE)
  stripchart(carcinomasTab, vertical = TRUE, las=2, pch = 20, col=c("olivedrab","olivedrab", "salmon","salmon", "royalblue","royalblue", "goldenrod","goldenrod"))
  boxplot(carcinomasTab, add=TRUE, xaxt='n', yaxt='n')
  text(x = 1, y=1, labels = sigAcomparison, cex = 0.5)
  text(x = 3, y=0.95, labels = sigBcomparison, cex = 0.5)
  text(x = 5, y=0.9, labels = sigCcomparison, cex = 0.5)
  text(x = 7, y=0.85, labels = sigDcomparison, cex = 0.5)
dev.off()


#stats for adenomas
adenomaTab <- plotTab[c(12:16),]

sigAcomparison <- wilcox.test(adenomaTab[["A"]], adenomaTab[["Ab"]])$p.value
sigBcomparison <- wilcox.test(adenomaTab[["B"]], adenomaTab[["Bb"]])$p.value
sigCcomparison <- wilcox.test(adenomaTab[["C"]], adenomaTab[["Cb"]])$p.value
sigDcomparison <- wilcox.test(adenomaTab[["D"]], adenomaTab[["Db"]])$p.value

pdfName <- paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-boxplots.adenomas.pdf", sep="")
pdf(file=pdfName, onefile=TRUE, width=5, height=5)
par(mar=c(7,5,5,5), xpd=TRUE)
stripchart(adenomaTab, ylim = c(0,0.8), vertical = TRUE, las=2, pch = 20, col=c("olivedrab","olivedrab", "salmon","salmon", "royalblue","royalblue", "goldenrod","goldenrod"))
boxplot(adenomaTab, add=TRUE, xaxt='n', yaxt='n')
text(x = 1, y=1, labels = sigAcomparison, cex = 0.5)
text(x = 3, y=0.95, labels = sigBcomparison, cex = 0.5)
text(x = 5, y=0.9, labels = sigCcomparison, cex = 0.5)
text(x = 7, y=0.85, labels = sigDcomparison, cex = 0.5)
dev.off()


#stats for Lynch and MSI
LynchTab <- plotTab[c(4, 17:19),]

sigAcomparison <- wilcox.test(LynchTab[["A"]], LynchTab[["Ab"]])$p.value
sigBcomparison <- wilcox.test(LynchTab[["B"]], LynchTab[["Bb"]])$p.value
sigCcomparison <- wilcox.test(LynchTab[["C"]], LynchTab[["Cb"]])$p.value
sigDcomparison <- wilcox.test(LynchTab[["D"]], LynchTab[["Db"]])$p.value

pdfName <- paste(sampleList[1,6], holdingDir, "162906.EMu/Assigned-boxplots.lynch.pdf", sep="")
pdf(file=pdfName, onefile=TRUE, width=5, height=5)
par(mar=c(7,5,5,5), xpd=TRUE)
stripchart(LynchTab, ylim = c(0,0.8), vertical = TRUE, las=2, pch = 20, col=c("olivedrab","olivedrab", "salmon","salmon", "royalblue","royalblue", "goldenrod","goldenrod"))
boxplot(LynchTab, add=TRUE, xaxt='n', yaxt='n')
text(x = 1, y=1, labels = sigAcomparison, cex = 0.5)
text(x = 3, y=0.95, labels = sigBcomparison, cex = 0.5)
text(x = 5, y=0.9, labels = sigCcomparison, cex = 0.5)
text(x = 7, y=0.85, labels = sigDcomparison, cex = 0.5)
dev.off()


################## 2. compare raw signatures ################





