# Differential Isoform Analysis
This repository details the pipeline used to analyse differential isoforms in fetal and adult cardiac fibroblasts. The aim of this project is to identify important isoforms and exons that contribute to the switch between fetal and adult cardiac fibroblasts.

The general project workflow is as follows: 
1. Transcript abundaces were quantified using raw FASTQ files of scRNA-seq data using Kallisto. 
2. Differential isoform analysis was performed on the Kallisto output using Sleuth. 
3. Differential exon usage of isoforms were analysed using DEXSeq. 
