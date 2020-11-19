library(sleuth)
library(data.table)
library(dplyr)
library(splitstackshape)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(gplots)
library(RColorBrewer)

##Sleuth Differential Isoform Analysis
#Set file directory
base_dir <- "/Users/a.h./Documents/GIS/Fibro Differential"
sample_id <- dir(file.path(base_dir, "Kallisto"))
kal_dirs <- file.path(base_dir, "Kallisto", sample_id)

#Load a table that describes experimental condition and path of each sample
condition <- c('fHCF', 'fHCF', 'fHCF', 'aHCF','aHCF','aHCF')
s2c <- data.table("sample" = sample_id, "condition" = condition, "path" = kal_dirs)

#Convert to sleuth object
sleuth_obj <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, ~condition)

#Fit full model 
full_model <- sleuth_fit(sleuth_obj, ~condition, 'full')
reduced_model <- sleuth_fit(full_model, ~1, 'reduced')

#Perform a second "reduced model" to identify differential expressed transcripts
likelihood_test <- sleuth_lrt(reduced_model, "reduced", "full")
#Check intercept name to perform wald_test
models(likelihood_test)
wald_test <- sleuth_wt(reduced_model, paste0('conditionfHCF'), 'full')

#Examine sleuth tests. Take note that WT's beta values are relative to aHCF.
results_lrt_table <- sleuth_results(likelihood_test, 'reduced:full', 'lrt')
#sleuth_significant <- dplyr::filter(results_table, qval <= 0.05)
results_wald_table <- sleuth_results(wald_test, 'conditionfHCF', test_type = 'wt')

#Retrieving avg tpm values for fHCF and aHCF for each of the differential isoform
RH072_tpm <- sleuth_obj$kal[[1]]$abundance
RH072_tpm <- RH072_tpm[match(sleuth_significant$target_id,RH072_tpm$target_id),]
RH073_tpm <- sleuth_obj$kal[[2]]$abundance
RH073_tpm <- RH073_tpm[match(sleuth_significant$target_id,RH073_tpm$target_id),]
RH074_tpm <- sleuth_obj$kal[[3]]$abundance
RH074_tpm <- RH074_tpm[match(sleuth_significant$target_id,RH074_tpm$target_id),]
fHCF_tpm <- data.frame(cbind(RH072 = RH072_tpm$tpm, RH073 = RH073_tpm$tpm, RH074 = RH074_tpm$tpm))
rownames(fHCF_tpm) <- RH072_tpm$target_id
fHCF_tpm$avg_tpm <- apply(fHCF_tpm, 1, function(x) log10(mean(x)+1))

RH075_tpm <- sleuth_obj$kal[[4]]$abundance
RH075_tpm <- RH075_tpm[match(sleuth_significant$target_id,RH075_tpm$target_id),]
RH076_tpm <- sleuth_obj$kal[[5]]$abundance
RH076_tpm <- RH076_tpm[match(sleuth_significant$target_id,RH076_tpm$target_id),]
RH077_tpm <- sleuth_obj$kal[[6]]$abundance
RH077_tpm <- RH077_tpm[match(sleuth_significant$target_id,RH077_tpm$target_id),]
aHCF_tpm <- data.frame(cbind(RH075 = RH075_tpm$tpm, RH076 = RH076_tpm$tpm, RH077 = RH077_tpm$tpm))
rownames(aHCF_tpm) <- RH075_tpm$target_id
aHCF_tpm$avg_tpm <- apply(aHCF_tpm, 1, function(x) log10(mean(x)+1))
sleuth_significant <- cbind(sleuth_significant, fHCF_avgTpm = fHCF_tpm$avg_tpm, aHCF_avgTpm = aHCF_tpm$avg_tpm)

#Calculate avg_logFC, avg_tpm and -log10(FDR) for each DI (Additional info for plots)
sleuth_significant$avg_logFC <- log2(sleuth_significant$aHCF_avgTpm +1) - log2(sleuth_significant$fHCF_avgTpm+1)
sleuth_significant$avg_tpm <- (sleuth_significant$aHCF_avgTpm + sleuth_significant$fHCF_avgTpm)/2
sleuth_significant$neglog_FDR <- -log10(sleuth_significant$qval)

#Assigning columns to each of the IDs
sleuth_significant <- cSplit(sleuth_significant, "target_id", "|")
sleuth_significant <- rename(sleuth_significant, target_id = target_id_1, gene_id = target_id_2, 
                             vega_gene_id = target_id_3, vega_transcript_id = target_id_4, 
                             transcript_name = target_id_5, gene_name = target_id_6, 
                             transcript_length = target_id_7, transcript_type = target_id_8)

#Plots/Matrices
conditions_scatter <- ggplot(sleuth_significant, aes(fHCF_avgTpm, aHCF_avgTpm)) + 
  geom_point() + theme(panel.background = element_blank(), axis.line = element_line()) + 
  geom_text_repel(size = 2, aes(label = ifelse((fHCF_avgTpm < 1 & aHCF_avgTpm > 2.1)|(fHCF_avgTpm < 2.5 & aHCF_avgTpm > 3)|(fHCF_avgTpm > 2.3 & aHCF_avgTpm < 1)|(fHCF_avgTpm > 3 & aHCF_avgTpm < 2.1), as.character(gene_name), '')), hjust = 0, vjust = 0.5, segment.colour = NA)

MA_plot <- ggplot(sleuth_significant, aes(avg_tpm, avg_logFC)) + geom_point() +
  theme(panel.background = element_blank(), axis.line = element_line()) +
  geom_text_repel(size = 2,aes(label = ifelse((avg_tpm > 0.6 & avg_logFC > 1.3)|(avg_tpm > 0.6 & avg_logFC < -1.3)|(avg_tpm > 1.2 & avg_logFC < -1)|(avg_tpm > 1.2 & avg_logFC > 1), as.character(gene_name), '')), hjust = 0, vjust = 1, segment.color = NA)

Volcano_plot <- ggplot(sleuth_significant, aes(avg_logFC, neglog_FDR)) + geom_point() + 
  theme(panel.background = element_blank(), axis.line = element_line()) +
  geom_text_repel(size = 2, aes(label = ifelse(neglog_FDR > 4.5, as.character(gene_name), '')),hjust = 1, vjust = -1, segment.colour = NA)

RH072_tpm$tpm <- log10(RH072_tpm$tpm +1)
RH073_tpm$tpm <- log10(RH073_tpm$tpm +1)
RH074_tpm$tpm <- log10(RH074_tpm$tpm +1)
RH075_tpm$tpm <- log10(RH075_tpm$tpm +1)
RH076_tpm$tpm <- log10(RH076_tpm$tpm +1)
RH077_tpm$tpm <- log10(RH077_tpm$tpm +1)
GeneTpm_Matrix <- data.frame(cbind(RH072 = RH072_tpm[1:30,]$tpm, RH073 = RH073_tpm[1:30,]$tpm, RH074 = RH074_tpm[1:30,]$tpm, RH075 = RH075_tpm[1:30,]$tpm, RH076 = RH076_tpm[1:30,]$tpm, RH077 = RH077_tpm[1:30,]$tpm))
GeneTpm_Matrix$Transcript_id <- sleuth_significant[1:30,]$transcript_id
GeneTpm_Matrix <- cbind(data.frame(Gene = sleuth_significant[1:30,]$gene_name), GeneTpm_Matrix)
GeneTpm_Matrix <- GeneTpm_Matrix[order(GeneTpm_Matrix$Gene),]

Heatmap_matrix <- aHCF_tpm[match(rownames(fHCF_tpm), rownames(aHCF_tpm)),]
Heatmap_matrix <- as.matrix(cbind(fHCF_tpm, aHCF_tpm))
correlation_matrix <- 1-cor(Heatmap_matrix)
clust <- hclust(as.dist(correlation_matrix), method = "ward.D2")
Heatmap_plot <- heatmap.2(as.matrix(correlation_matrix), trace = 'none', density.info = 'none', Rowv = as.dendrogram(clust), labRow = FALSE, dendrogram = 'none')

pathway_top500_DI <- sleuth_significant[1:500]$gene_name
pathway_all_DI <- sleuth_significant$gene_name

top30_genes_isoform <- unique(sleuth_significant$gene_name)[1:30]
isoform_FC  <- sleuth_significant[sleuth_significant$gene_name %in% top30_genes_isoform,]
count_isoforms <- isoform_FC %>% group_by(gene_name) %>% summarise(n_distinct(transcript_id))
count_isoforms <- count_isoforms[rev(order(count_isoforms$`n_distinct(transcript_id)`)),]
gene_isoformFC_data <- data.table(Gene = isoform_FC$gene_name, Transcript = isoform_FC$transcript_name, avg_logFC = isoform_FC$avg_logFC)
gene_isoformFC_data <- gene_isoformFC_data[order(match(gene_isoformFC_data$Gene, count_isoforms$gene_name)),]
gene_isoformFC_data$Gene <- factor(gene_isoformFC_data$Gene, levels = rev(count_isoforms$gene_name))
genetoisoform_avglogFC <- ggplot(gene_isoformFC_data, aes(avg_logFC, Gene)) + geom_point() + 
  theme(panel.background = element_blank(), axis.line = element_line(), axis.ticks.y = element_blank())

DG_paper <- fread("/Users/a.h./Documents/GIS/Fibro Differential/Differential_Genes_Paper.diff")
significant_DGs <- DG_paper[DG_paper$q_value <= 0.05,]
new_DGs <- unique(setdiff(pathway_all_DI, significant_DGs$gene))
new_DGs_info <- sleuth_significant[sleuth_significant$gene_name %in% new_DGs,]
new_DGs_info <- new_DGs_info[rev(order(new_DGs_info$qval)),]
top30_newDGs_logFC <- new_DGs_info$gene_name[1:30]
top30_newDGs_info <- new_DGs_info[new_DGs_info$gene_name %in% top30_newDGs_logFC,]
top30_newDGs_info <- top30_newDGs_info[1:30,]
top30_newDGs_info <- data.frame(Gene = top30_newDGs_info$gene_name, Isoform = top30_newDGs_info$transcript_id, avg_logFC = top30_newDGs_info$avg_logFC)

#####Sleuth P-Value Aggregation Analysis
#Set file directory
base_dir <- "/Users/a.h./Documents/GIS/Fibro Differential"
sample_id <- dir(file.path(base_dir, "Kallisto_Aggregate"))
kal_dirs <- file.path(base_dir, "Kallisto_Aggregate", sample_id)

#Load a table that describes experimental condition and path of each sample
condition <- c('fHCF', 'fHCF', 'fHCF', 'aHCF','aHCF','aHCF')
sc <- data.table("sample" = sample_id, "condition" = condition, "path" = kal_dirs)

#Set a target to gene table
kallisto_result <- fread("/Users/a.h./Documents/GIS/Fibro Differential/Kallisto/RHf072/abundance.tsv")
id_column <- kallisto_result[,1]
id_column <- cSplit(id_column, "target_id", "|")
id_column <- rename(id_column, target_id = target_id_1, gene_id = target_id_2, 
                             vega_gene_id = target_id_3, vega_transcript_id = target_id_4, 
                             transcript_name = target_id_5, gene_name = target_id_6, 
                             transcript_length = target_id_7, transcript_type = target_id_8)
id_column <- id_column[,c("target_id", "gene_id", "gene_name")]

#Set up sleuth obj
sleuth_pval_obj <- sleuth_prep(sc, extra_bootstrap_summary = TRUE, target_mapping = id_column, aggregation_column = 'gene_id', gene_mode = TRUE)

#Fit full model
model <- sleuth_fit(sleuth_pval_obj, ~condition ,'reduced')
model <- sleuth_fit(model, ~condition ,'full')

#Conduct differential gene analysis
lrt_model <- sleuth_lrt(model, 'reduced', 'full')
models(lrt_model)
wt_model <- sleuth_wt(model, paste0('conditionfHCF'), 'full')
  
#Conduct likelihood test for differential gene analysis
results_lrt <- sleuth_results(lrt_model, 'reduced:full', 'lrt', show_all = TRUE)
results_wald <- sleuth_results(wt_model, 'conditionfHCF', test_type = 'wt', show_all = TRUE)

#Gene-level TPM values
results <- sleuth_to_matrix(sleuth_pval_obj, "obs_norm", "tpm")


