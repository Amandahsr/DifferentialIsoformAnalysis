#This script details the pipeline for differential exon analysis using DEXSeq.
library(DEXSeq)

#Initialise parameters.
sample1_name <- ""
sample2_name <- ""
sample3_name <- ""
sample4_name <- ""
sample5_name <- ""
sample6_name <- ""
sample1_condition <- ""
sample2_condition <- ""
sample3_condition <- ""
sample4_condition <- ""
sample5_condition <- ""
sample6_condition <- ""

#Set file directory to exon counts files and gene annotation gff files. 
base_dir <- "counts folder"
countFiles <- list.files(base_dir, pattern = ".dexseq_count.txt$", full.names = TRUE)
flattenedFile <- list.files(base_dir, pattern = "gff$", full.names = TRUE)

#Create a table containing experimental condition of each sample.
condition <- c(sample1_condition, sample2_condition, sample3_condition, sample4_condition, sample5_condition, sample6_condition)
names <- c(sample1_name, sample2_name, sample3_name, sample4_name,sample5_name, sample6_name)
sampleTable <- data.frame(condition = condition)
rownames(sampleTable) <- names

#Initialise table as DEXSeq object.
dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData = sampleTable, 
                              design = ~ sample + exon + condition:exon, 
                              flattenedfile = flattenedFile)

#Normalise and estimate variability of data.
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)

#Test for model fit and estmate fold changes.
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

#Retrieve DEXSeq results.
dxr1 <- DEXSeqResults(dxd, independentFiltering = FALSE)

#Plot DEXSeq results of a particular gene.
plotDEXSeq(dxr1, "gene_id", legend = TRUE, expression = FALSE, 
           splicing = TRUE, cex.axis=0.8,cex=0.6, lwd=2)

