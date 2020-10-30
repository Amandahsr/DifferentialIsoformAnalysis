# Differential Isoform Analysis
This repository details the pipeline used to analyse differential isoforms in fetal and adult cardiac fibroblasts. The aim of this project is to identify important isoforms and exons that contribute to the switch between fetal to adult cardiac fibroblasts.

The general project workflow is as follows: Transcript abundaces were quantified using raw FASTQ files of scRNA-seq data using Kallisto. Differential isoform analysis was performed on the Kallisto output using Sleuth. Lastly, differential exon usage of isoforms were analysed using DEXSeq. The contribution of splicing factors to isoform switches was were also analysed by identifying splicing factor binding sites using SFMap and SpliceAID. 
