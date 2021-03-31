### Passed Global Arguments ###

args <- c("/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/scripts", "/pellmanlab/stam_niko/data/processed_bam")
scriptsdir <- args[1]
dirpath <- args[2]

### Data Read-ins ###

controlSampleIDs2 <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs2.rds", wkpath))
controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", wkpath))
#geneSpecific <- readRDS(file = sprintf("%s/aggregated_results/geneSpecific.rds", wkpath))
rsemtpm <- readRDS(file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", wkpath))
ASE <- readRDS(file = sprintf("%s/aggregated_results/ASE.rds", wkpath))
#coding <- readRDS(sprintf("%s/aggregated_results/coding_snps.rds", wkpath))
#configs <- readRDS(sprintf("%s/aggregated_results/param_config_list.rds", wkpath))
geneRanges <- readRDS(sprintf("%s/aggregated_results/geneRanges.rds", wkpath))
#adt.old <- readRDS("/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/adt.na.rds")
alleles <- readRDS(sprintf("%s/aggregated_results/alleles.all.rds", wkpath))
alleles.frac <- readRDS(sprintf("%s/aggregated_results/alleles.all.frac.rds", wkpath))

############################################################################
# Data Restructuring of Variants and Gene Expr (TPM) into format for Models:
############################################################################

#Prepare data for GSSM: features totalExp, A, B, minor allele freq, gene interdist
#------------------------------------------------------------------------------------

#data preparation for CNV estimates:
#----------------------------------
source('/homes10/ejacob/WORK/secondaryanalysis/methods_paper_results/R/CNML.R')

#Total expression with +1:
#-------------------------
source('/pellmanlab/nikos/Stam_Etai_Scripts/fromRawTPMtoExprsRatio.R')
source('/pellmanlab/nikos/Stam_Etai_Scripts/generate_baseline_datasets_for_ML_and_stats_analysis.R')
tmppath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2"
configs <- readRDS(sprintf("%s/data/param_config_list.rds", tmppath))

generate_adt_adt0_adt.na_and_nonzeros_data()

