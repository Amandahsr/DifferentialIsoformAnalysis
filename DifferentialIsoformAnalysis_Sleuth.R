#This script details the pipeline for differential isoform analysis using Sleuth.
#Bootstrapping must be performed in Kallisto for Sleuth analysis. 

library(sleuth)
library(data.table)

##1. Process isoform data.
##########################################################
#Set file directory to Kallisto results. 
#Kallisto results for each sample should be kept in its own folder.
#Folder structure: "Kallisto" folder -> Sample folders -> Kallisto results (each sample folder)
base_dir <- "Kallisto Results"
sample_id <- dir(file.path(base_dir, "Kallisto"))
kal_dirs <- file.path(base_dir, "Kallisto", sample_id)

#Create a table containing experimental condition and file directory of each sample.
condition <- c('sample1_condition', 'sample2_condition', 'sample3_condition', 'sample4_condition','sample5_condition','sample6_condition')
s2c <- data.table("sample" = sample_id, "condition" = condition, "path" = kal_dirs)

#Initialise table as sleuth object.
sleuth_obj <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, ~condition)

#Fit results into a model for differential analysis.
full_model <- sleuth_fit(sleuth_obj, ~condition, 'full')
reduced_model <- sleuth_fit(full_model, ~1, 'reduced')

##2. Peform differential isoform analysis. Choose either likelihood test or wald test.
##########################################################
#2a. Likelihood test do not provide estimated fold change.
likelihood_test <- sleuth_lrt(reduced_model, "reduced", "full")

#Fold change can be estimated as avglogFC = log2(avgTpm("Condition")) - log2(avgTpm("Control")). 
fHCF_Tpm <- data.frame(cbind(sample1 = sample1_Tpm, sample2 = sample2_Tpm, sample3 = sample3_Tpm))
fHCF_avgTpm <- apply(fHCF_Tpm, 1, function(x) mean(x))
aHCF_Tpm <- data.frame(cbind(sample4 = sample4_Tpm, sample5 = sample5_Tpm, sample6 = sample6_Tpm))
aHCF_avgTpm <- apply(aHCF_Tpm, 1, function(x) mean(x))  
avglogFC <- log2(aHCF_avgTpm +1) - log2(fHCF_avgTpm+1)

#2b. Wald test provides estimated fold change (B value). 
models(reduced_model) #Extract the model's intercept.
wald_test <- sleuth_wt(reduced_model, paste0('intercept'), 'full')

##3. Extract differential analysis results.
##########################################################
#3a. Likelihood test results.
lrt_results <- sleuth_results(likelihood_test, 'reduced:full', 'lrt')

#3b. Wald test results.
wald_results <- sleuth_results(wald_test, 'intercept', test_type = 'wt')
