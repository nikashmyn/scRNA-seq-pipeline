
####### Global Variables ########
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025misc", "SIS1025targ")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
wkpath <- dirpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", wkpath)
system(mk_agg_dir)

require(data.table)
require(readxl)
source( sprintf('%s/scripts/runRNAseqPipelineFor_BF9VP.utils.R', scriptsdir) )
source( sprintf('%s/scripts/generate_baseline_datasets_for_ML_and_stats_analysis.R', scriptsdir) )
source( sprintf('%s/plots/class_prob_to_cn_level.R', scriptsdir) )

#new annotation lists are not comprehensive. So older versions are needed for all samples.
#anno <- data.table(readRDS( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))
#anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_11_6_21.csv", datadir)))
anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir)))

#centromeres object is called inside AddFeaturesandPhasing function
centromeres <- readRDS( sprintf("%s/centromeres.rds", datadir) )

#######################################################################
#Collection of total expression data (controls and samples to predict):
#######################################################################

# Account for the special case where SISf_L1 and L2 are the same cells
SISf_L1 <- "SIS1025f_Lane1"
SISf_L2 <- "SIS1025f_Lane2"
#This will be changed to 1 if SISf L1 and L2 are being called
SISf_Switch <- 0

if (SISf_L1 %in% args) {
  if (SISf_L2 %in% args) {
    #Average the expression and make new object
    L1_tpm <- readRDS(sprintf("%s/%s/expr_results/%s_rsemtpm.rds", wkpath, SISf_L1, SISf_L1))
    L2_tpm <- readRDS(sprintf("%s/%s/expr_results/%s_rsemtpm.rds", wkpath, SISf_L2, SISf_L2))
    SISf_rsemtpm <- (L1_tpm + L2_tpm)/2
    
    #Add the variants and make new object
    L1_var <- readRDS(sprintf("%s/%s/variant_results/%s_raw_variants.rds", wkpath, SISf_L1, SISf_L1))
    L2_var <- readRDS(sprintf("%s/%s/variant_results/%s_raw_variants.rds", wkpath, SISf_L2, SISf_L2))
    SISf_variants <- list(A = L1_var$A + L2_var$A, B = L1_var$B + L2_var$B, features = L1_var$features)
   
    #get QC information for SISf lanes
    L1_QC <- read.csv(sprintf("%s/%s/QC/%s_QC.csv", wkpath, SISf_L1, SISf_L1))
    L2_QC <- read.csv(sprintf("%s/%s/QC/%s_QC.csv", wkpath, SISf_L2, SISf_L2))
    SISf_QC <- cbind( L1_QC[,c(1)], round((L1_QC[,-c(1)] + L2_QC[,-c(1)])/2) ) #average the QC between the two lanes
    colnames(SISf_QC) <- colnames(L1_QC) #same col names this can be rbinded later
    
    #Remove the paths from the rsem_paths object
    experiments <- experiments[experiments != SISf_L1]
    experiments <- experiments[experiments != SISf_L2]
    
    #Flip the SISf switch 
    SISf_Switch <- 1
  }
}

#Get paths to all TPM objects
rsem_paths <- list()
for (i in experiments) {rsem_paths <- cbind(rsem_paths, sprintf("%s/%s/expr_results/%s_rsemtpm.rds", wkpath, i, i))}

  
# Now we will bind (column-wise) all the experiments together 
all_rsemtpm = list()
for (j in rsem_paths) {
  if (length(all_rsemtpm) != 0) {
    tmp <- readRDS(j)
    stopifnot(sum(rownames(tmp) != rownames(all_rsemtpm)) == 0)
    all_rsemtpm <- cbind(all_rsemtpm, tmp)
  } else {
    all_rsemtpm <- readRDS(j)
  }
}

# Now bind the special case if switch in on
if (SISf_Switch == 1) {
  all_rsemtpm <- cbind(all_rsemtpm, SISf_rsemtpm)
}

tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(i)  {
  if (!is.na(anno$Fastq_files[i])) {
    res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[i]), ignore.case = T, value=T, x = colnames(all_rsemtpm))
    if (length(res) < 1) {  res <- colnames(all_rsemtpm)[which(colnames(all_rsemtpm) %in% anno$Fastq_files[i])] }
    data.table(raw_ID = res, ID = anno$Fastq_files[i], name = anno$WTA.plate[i])}}))
all_rsemtpm_ids <- data.table(t(all_rsemtpm), keep.rownames = "raw_ID")
setkeyv(tmp, "raw_ID")
setkeyv(all_rsemtpm_ids, "raw_ID")
common_samples <- intersect(tmp$raw_ID, all_rsemtpm_ids$raw_ID)
rsemtpm <- as.data.frame(cbind(tmp[which(tmp$raw_ID %in% common_samples),], all_rsemtpm_ids[which(all_rsemtpm_ids$raw_ID %in% common_samples),]))
rownames(rsemtpm) <- rsemtpm$name
rsemtpm <- t(rsemtpm[,-c(1:4)])
rsemtpm <- as.data.frame(rsemtpm)
message(sprintf("The dimensions of the TPM object is %s %s \n", dim(rsemtpm)[1], dim(rsemtpm)[2]))

saveRDS(object = rsemtpm, file = sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", wkpath))
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", wkpath))


############################################
# Control Samples Collection for Annotations
############################################

controlSampleIDs <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs <- controlSampleIDs[controlSampleIDs %in% colnames(rsemtpm)]
saveRDS(controlSampleIDs, sprintf("%s/aggregated_results/controlSampleIDs.rds", wkpath))

controlSampleIDs2 <- anno[ (mainGroup == "A_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs2 <- controlSampleIDs2[controlSampleIDs2 %in% colnames(rsemtpm)]
saveRDS(controlSampleIDs2, sprintf("%s/aggregated_results/controlSampleIDs2.rds", wkpath))

#############################################################################################################################
#Collection of Variant data (phased/biallelic)
#############################################################################################################################

#Merge current run with previous data:
#-------------------------------------
#loading previous allele data (already phased):
#These data also includes biallelicFilter based on the function: loadAndPhaseSNPdataFromRPE1ExperimentsBasedOnHaplotypeAndExpression.R in workflows
variants_paths <- list()
for (i in experiments) {variants_paths <- cbind(variants_paths, sprintf("%s/%s/variant_results/%s_raw_variants.rds", wkpath, i, i))}

all_vars = list()
for (j in variants_paths) {
  if (length(all_vars) != 0) {
    tmp <- readRDS(j)
    stopifnot(sum(rownames(tmp$features) != rownames(all_vars$features)) == 0)
    all_vars$A <- cbind(all_vars$A, tmp$A)
    all_vars$B <- cbind(all_vars$B, tmp$B)
  } else {
    all_vars <- readRDS(j)
  }
}

if (SISf_Switch == 1) {
  all_vars$A <- cbind(all_vars$A, SISf_variants$A)
  all_vars$B <- cbind(all_vars$B, SISf_variants$B)
  
}


saveRDS(object = all_vars, file = sprintf("%s/aggregated_results/all_vars.rds", wkpath))
all_vars <- readRDS(sprintf("%s/aggregated_results/all_vars.rds", wkpath))

#reduces mym sample names from filename to anno list name. Currently checking all anno lists
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(i)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[i]), ignore.case = T,  value=T, x = colnames(all_vars$A))
  if (length(res) < 1) {  res <- colnames(all_vars$A)[which(colnames(all_vars$A) %in% anno$Fastq_files[i])] }
  data.table(raw_ID = res, ID = anno$Fastq_files[i], name = anno$WTA.plate[i])}))
all_vars_ids <- colnames(all_vars$A)
setkeyv(tmp, "raw_ID")
common_samples <- intersect(tmp$raw_ID, all_vars_ids)
common_IDs <- tmp$name[which(tmp$raw_ID %in% common_samples)]
common_IDs <- as.character(common_IDs)
all_vars$A <- all_vars$A[,..common_samples]
setnames(all_vars$A, common_IDs)
all_vars$B <- all_vars$B[,..common_samples]
setnames(all_vars$B, common_IDs)

#phase SNPs:
# there are duplicate names from double renaming. possibly phase pre merge. 
source( sprintf('%s/scripts/phaseMyAllelesBasedOnHaplotypes.R', scriptsdir) )
mym <- phaseMyAllelesBasedOnHaplotypes(m = all_vars)
saveRDS(object = mym, file = sprintf("%s/aggregated_results/haplotypePhasedADs.rds", wkpath))
mym <- readRDS(file = sprintf("%s/aggregated_results/haplotypePhasedADs.rds", wkpath))

#Add features to snps
ASE <- AddFeaturesandPhasing(mym, anno = anno, centromeres = centromeres)
saveRDS(object = ASE, file = sprintf("%s/aggregated_results/ASE.rds", wkpath))
ASE <- readRDS(file = sprintf("%s/aggregated_results/ASE.rds", wkpath))

#Allele level data:
alleles.all <- prepareAlleleData(m = ASE, useFraction = F , controlSampleIDs = controlSampleIDs)
saveRDS(alleles.all, file = sprintf("%s/aggregated_results/alleles.all.rds", wkpath))

alleles.all.frac <- prepareAlleleData(m = ASE, useFraction = T , controlSampleIDs = controlSampleIDs)
saveRDS(alleles.all.frac, file = sprintf("%s/aggregated_results/alleles.all.frac.rds", wkpath))

###################################
# Aggregate QC for all experiments:
###################################
QC_paths <- list()
for (i in experiments) {QC_paths <- cbind(QC_paths, sprintf("%s/%s/QC/%s_QC.csv", wkpath, i, i))}

all_QC = list()
for (j in QC_paths) {
  if (length(all_QC) != 0) {
    tmp <- read.csv(j)
    all_QC <- rbind(all_QC, tmp)
  } else {
    all_QC <- read.csv(j)
  }
}

if (SISf_Switch == 1) {
  all_QC <- rbind(all_QC, SISf_QC)
}
colnames(all_QC) <- c("id","N_unmapped","N_multimapping","N_noFeature","N_ambiguous","geneCounts","th1","th5","th10","totalReads")
all_QC <- all_QC[which(all_QC$id %in% anno$WTA.plate),]  #Exclude cells not in the annotation list
saveRDS(object = all_QC, file = sprintf("%s/aggregated_results/all_QC.rds", wkpath))
write.csv(all_QC, file = sprintf("%s/aggregated_results/all_QC.csv", wkpath))
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", wkpath))

hist(all_QC$th5, xlab = "Genes with >5 read coverage", main = "Th5 Distribution Across All Samples")
abline(v = 4000)
cutoff_vals <- quantile(all_QC$th5, c(.05, .10, .25))

#######################
# Organize Data for ML:
#######################

#Standard set of annotations for all genes
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))

#list of preset parameters for the experiment
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))
configs$minDetectionLevel <- 50

#To me it looks like etai's df is standardized as well
# Use the first for visuals and second for ML. Curretly ran second version.
#generate_adt_adt0_adt.na_and_nonzeros_data(dirpath = wkpath, th = 4, scriptsdir = scriptsdir, normBySd = T) 
generate_adt_adt0_adt.na_and_nonzeros_data(dirpath = wkpath, th = 4, scriptsdir = scriptsdir, datadir = datadir, normBySd = F) 

#creates coding_snp rds object for ML
generate_coding_snps_data(dirpath = wkpath)

#creates ASE.coding.rds and ASE.noncoding.rds
get_snps_fraction_bins_zscore(saveMe = T, dirpath = wkpath)

#creates the nonzeros.zs.bin50.rds and nonzeros.bin50.rds objects
nonzeros <- compute_nonZero_bins_zscore(saveMe = T, dirpath = wkpath)
#This could be reexported if needed

#Creates arm level annotations
generate_ganno2(dirpath = wkpath, datapath = datadir)
#genespecific could be created as well. not needed for now. function needs to be fixed.
#generate_RPE1_GeneSpecificData(dirpath = wkpath)

#Create chr binned normalized matrices for TPM and Variant Counts
source(sprintf("%s/scripts/normalize_tpm_and_vars_bychr.R", scriptsdir))

#Generate the arms object made of SSM, MA, & allelic CN predictions
datapath <- datadir
CNdir <- sprintf("%s/CN_data", dirpath)
mk_CN_dir <- sprintf("mkdir -p %s/CN_data", dirpath)
system(mk_CN_dir)
source(sprintf("%s/scripts/run_all_CN_predictions_for_sample_nikos.R", scriptsdir))

#####################################
# Add QC to annotations for new list:
#####################################

#Reduce to mergable dataframes
anno_QC <- data.table(anno[which(anno$WTA.plate %in% all_QC$id),c(1:22)])
QC_anno <- data.table(all_QC[which(all_QC$id %in% anno$WTA.plate),])

#Get the sum of all variants for each cell by allele
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))
tmp_QC_coding_A <- data.table(cbind(WTA.plate = colnames(coding$cnts.A[,-c(1:2)]), var_counts_A = unlist(colSums(coding$cnts.A[,-c(1:2)]))))
tmp_QC_coding_B <- data.table(cbind(WTA.plate = colnames(coding$cnts.B[,-c(1:2)]), var_counts_B = unlist(colSums(coding$cnts.B[,-c(1:2)]))))
QC_coding_A <- tmp_QC_coding_A[which(tmp_QC_coding_A$WTA.plate %in% QC_anno$id),]
QC_coding_B <- tmp_QC_coding_B[which(tmp_QC_coding_B$WTA.plate %in% QC_anno$id),]

#merge the dataframes
setkeyv(anno_QC, "WTA.plate")
setnames(QC_anno, old = "id", new = "WTA.plate")
setkeyv(QC_anno, "WTA.plate")
anno_w_QC <- merge(anno_QC, QC_anno, by="WTA.plate")
setkeyv(QC_coding_A, "WTA.plate")
setkeyv(QC_coding_B, "WTA.plate")
anno_w_QC_w_A <- merge(anno_w_QC, QC_coding_A, by = "WTA.plate")
new_anno <- merge(anno_w_QC_w_A, QC_coding_B, by = "WTA.plate")
write.csv(new_anno, sprintf("%s/aggregated_results/analysis_list.csv", dirpath))
#new_anno <- read.csv(sprintf("%s/aggregated_results/analysis_list.csv", dirpath))
#test <- new_anno[which(new_anno$var_counts_A < 1000 ),]
#quantile(new_anno$var_counts_A, .1)

################################################################################

#ready to run "Run_ML.R" script
print("Done")

