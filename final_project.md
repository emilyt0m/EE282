## Final Project Report

#### Emily Tom

#### EE282 Fall 2022

------------------------------------------------------------------------

### Introduction

The overall goal of this research project is to identify changes in the transcriptome in retinal pigment epithelial (RPE) cells during aging. To accomplish this, I will perform RNA sequencing on ARPE-19 cells, a widely used immortalized human RPE cell line, lacking *ELOVL2* (Elongation of Very Long Chain Fatty Acids Like-2). *Note: In my topic and analysis proposals, I had originally planned to perform RNA sequencing on RPE derived from human induced pluripotent stem cells (iPSCs) lacking* ELOVL2; *however, due to a delay from the Genomics Research and Technology Hub, I have decided to analyze this previously generated data set. However, the overall goal can still be accomplished, just in a different cell line.* A total of 6 samples (3 biological replicates in each group) will be analyzed.

For further background information on *ELOVL2* and my research project, please refer to my topic and analysis proposals.

------------------------------------------------------------------------

### Methods

#### Data Acquisition

Raw data was transferred from a previous lab member's data set to `data/raw` with the following command:

`cp -r /dfs4/som/qianlanx/RAW/RNA-Seq_siE2-ARPE19_04-2021 /pub/etom2/classrepos/EE282/data/raw`

#### Data Pre-Processing

**Read Quality**

FASTQC was used to analyze read quality. The script `fastqc.sub` located in the directory `/code/scripts` was used to run FASTQC on all `*.txt` FastQ files.

**Read Trimming**

Trimmomatic was used to remove adapter sequences and sequences with poor quality. The script `trim.sub` located in the directory `/code/scripts` was used to remove Illumina Universal adapters.

#### Reference Genome

The *Homo sapiens* GRCh38 genome was downloaded from HISAT2 using `wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz` and unzipped using `tar -xf grch38_genome.tar.gz`. This was a prebuilt index; however, if necessary, HISAT2 indexes can be built using `hisat2-build` on the desired reference fasta and gtf file.

#### Data Alignment

Data was aligned using HISAT2 using the script `hisat2.sub` located in the directory `/code/scripts`. Alignment summary files were saved in the directory `/alignment/alignment_summary`. The output `.sam` files were then converted to `.bam` files with `samtools`: `samtools view -bS *.sam > *.bam`.

#### Read Counting

`featureCounts` from the Subread package was used to count how many of the aligned reads map to genomic features, such as genes, exons, etc. The script `featurecount.sub` located in the directory `/code/scripts` was used and the reference gtf file was downloaded from HISAT2 using `wget https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz` and `tar -xf hg38_genome.tar.gz`.

#### Data Visualization and PCA Plotting

I used DESeq2 to create a sample-wise PCA plot by regularized log transformation (`rlog`) and variance stabilizing transformation (`varianceStabilizingTransformation`). The `rlog` transformation produces a similar variance stabilizing effect as variance stabilizing transformation, but `rlog` is more robust when size factors vary widely. This allows us to check for potential outliers that may need to be removed before proceeding with the following steps.

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

**Differentially Expressed Gene (DEG) Analysis**

DEG analysis was conducted with `DESeq2`, and the results were sorted by the adjusted p-value in ascending order. Results were subsetted for genes with Benjamini-Hochberg adjusted p-value\<0.05 and exported as `degresults.csv`.

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

**Gene Clustering and Heatmap**

Hierarchical clustering was performed for the top 25 variable genes using `heatmap.2` .

    #Gene clusters and heatmap
    library("RColorBrewer")
    library("gplots")
    ramp <- 1:3/3
    cols <- c(rgb(ramp,0,0), rgb(0,ramp,0), rgb(0,0,ramp), rgb(ramp,0,ramp))
    topVarGenes <- head(order(genefilter::rowVars(assay(rld)), decreasing=T), 25)
    heatmap.2(assay(rld)[topVarGenes,], scale="row", trace="none", density.info="none", dendrogram="none", margins=c(6,6), colv="FALSE", col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

**Volcano Plot**

A volcano plot was generated using `EnhancedVolcano` with a p-value cutoff = 10e-6.

    #Volcano Plot
    library(EnhancedVolcano)
    EnhancedVolcano(res, lab= rownames(res), x = 'log2FoldChange', y = 'pvalue')

**Functional Enrichment Analysis**

The DEGs were uploaded to the **D**atabase for **A**nnotation, **V**isualization and **I**ntegrated **D**iscovery ([`DAVID`](https://david.ncifcrf.gov/home.jsp)) for functional enrichment. Data was then uploaded to [SRplot](http://www.bioinformatics.com.cn/srplot) to generate figures.

------------------------------------------------------------------------

### Results

#### Read Quality

The FASTQC Reports showed good quality reads, with an average Phred score of 36.

#### Alignment Statistics

On average, 98.7% of the reads were aligned.

|     Sample      | \% Reads Aligned |
|:---------------:|:----------------:|
| nR136-L1-G2-P51 |       98.7       |
| nR136-L1-G2-P52 |      98.83       |
| nR136-L1-G2-P53 |      98.85       |
| nR136-L1-G2-P54 |      98.83       |
| nR136-L1-G2-P55 |      98.89       |
| nR136-L1-G2-P56 |      98.87       |

#### PCA Plot

**PCA Plot (rlog)**

![](output/figures/PCA_rld.png)

**PCA Plot (vsd)**

![](output/figures/PCA_vsd.png)

Both transformations produced similar PCA plots.

#### Gene Clustering and Heatmap

![](output/figures/heatmap.png)

**Volcano Plot**

![](output/figures/Volcano_Plot.png)

**Functional Enrichment Analysis**

## ![](output/figures/GO_CC.png)

------------------------------------------------------------------------

### Discussion

In this project, I have analyzed changes in the transcriptome in ARPE-19 cells lacking *ELOVL2*. I processed RNAseq data and generated plots to visualize the changes between the samples.

The PCA plots show that one CTRL and one KD sample were different from the other two samples in their respective groups. This is confirmed by the heatmap representation, as hierarchical clustering also grouped CTRL 1 and KD 1 further from the other two samples in their groups.

As expected, *ELOVL2* showed the most significant downregulation. Another significantly downregulated gene was *EIF3H*, or Eukaryotic translation initiation factor 3 subunit H, which contributes to translation initiation factor activity. This could then lead to the repression of synthesis of many other proteins.

*DOX1*, or Alpha-dioxygenase 1, was significantly upregulated in KD cells. *DOX1* is involved in protection against oxidative stress and cell death, which could be a protective mechanism by the cells. It is also responsible for catalyzing the primary oxygenation step of a variety of 14-20 carbon fatty acids, which could be compensating for the loss of *ELOVL2*, which elongates 20 and 22 carbon fatty acids.

GO analysis showed many cellular membrane- related changes in KD cells. I hypothesize that the products of ELOVL2, such as docosahexanoic acid (DHA) and other very long-chain fatty acids (VLC-PUFAs), are important for maintaining membrane fluidity and integrity, and therefore, loss of *ELOVL2* greatly disturbs the membrane composition of many cellular organelles.

Future directions of this project include analyzing RNAseq data from iPSC-RPE cells lacking *ELOVL2*, which was my initial goal. However, this preliminary analysis provided valuable experience in designing pipelines for my future projects.
