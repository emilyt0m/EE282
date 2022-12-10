#Install DESeq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

#Load counts table
countsTable <- read.csv("counts.csv", header=T, stringsAsFactors=FALSE)

rownames(countsTable) <- countsTable[,1]
countsTable <- countsTable[,-1]

library(DESeq2)
conds <- factor(c("CTRL","CTRL","CTRL","KD","KD","KD"))
colData = data.frame(condition = conds)
dds <- DESeqDataSetFromMatrix(countsTable, colData, design=~ condition)

# Exploratory data analysis
# Relationship between variance (sigma^2) and mean (mu) of genes across samples
mu = rowMeans(countsTable)
sigma2 = apply(countsTable,1,var)
plot(log(mu), log(sigma2), xlim = c(0, log(max(mu))), ylim = c(0, log(max(mu))), pch=16, cex=0.3)

abline(0, 1, col="red")

#Correct for heteroskedasticity
rld <- rlog(dds)
colnames(rld) = colnames(countsTable)

plotPCA(rld, intgroup = "condition")

vsd <- varianceStabilizingTransformation(dds)
colnames(vsd) = colnames(countsTable)

plotPCA(vsd, intgroup = "condition")

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
#Meaning of columns
mcols(res, use.names=TRUE)

summary(res)

#Sort results by adjusted p value in ascending order
res <- res[order(res$padj),]

#Subset results for genes with BH adjusted p value < 0.05
deg = subset(res, padj < 0.05)

write.csv(as.data.frame(res), file = "degresults.csv")

#Gene clusters and heatmap
library("RColorBrewer")
library("gplots")
ramp <- 1:3/3
cols <- c(rgb(ramp,0,0), rgb(0,ramp,0), rgb(0,0,ramp), rgb(ramp,0,ramp))
topVarGenes <- head(order(genefilter::rowVars(assay(rld)), decreasing=T), 25)
heatmap.2(assay(rld)[topVarGenes,], scale="row", trace="none", density.info="none", dendrogram="none", margins=c(6,6), colv="FALSE", col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

#Volcano Plot
library(EnhancedVolcano)
EnhancedVolcano(res, lab= rownames(res), x = 'log2FoldChange', y = 'pvalue')
