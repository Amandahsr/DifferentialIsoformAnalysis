library(DEXSeq)

#Setup file paths
base_dir <- "/Users/a.h./Documents/GIS/Fibro Differential/dexseq_counts"
countFiles <- list.files(base_dir, pattern = ".dexseq_count.txt$", full.names = TRUE)
flattenedFile <- list.files(base_dir, pattern = "gff$", full.names = TRUE)

#Setup sample table
condition <- c('fHCF', 'fHCF', 'fHCF', 'aHCF','aHCF','aHCF')
names <- c("RHF072", "RHF073","RHF074","RHF075","RHF076","RHF077")
sampleTable <- data.frame(condition = condition)
rownames(sampleTable) <- names

#Prepare DEXSeq object
dxd <- DEXSeqDataSetFromHTSeq(countFiles, sampleData = sampleTable, 
                              design = ~ sample + exon + condition:exon, 
                              flattenedfile = flattenedFile)

#Normalise, and estimate variability of data
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)

#Test for model fit, and estmate fold changes
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

#Retrieve FN1 results
dxr1 <- DEXSeqResults(dxd, independentFiltering = FALSE)
plotDEXSeq(dxr1, "ENSG00000115414.14", legend = TRUE,expression = FALSE, splicing = TRUE, cex.axis=0.8,cex=0.6, lwd=2)

