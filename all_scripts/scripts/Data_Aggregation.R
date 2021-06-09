
####### Global Variables ########
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b",  "SIS1025d", "SIS1025e", "SIS1025f_Lane1" , "SIS1025f_Lane2")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
wkpath <- args[2]
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
anno <- data.table(readRDS( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))

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

#These change the names of the columns from their file names to the sample name. 
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[x]), ignore.case = T, x = colnames(all_rsemtpm))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
colnames(all_rsemtpm)[tmp[idx>0]$idx] <- tmp[idx>0]$id
print(all_rsemtpm[c(1),c(1)])
saveRDS(all_rsemtpm, sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", wkpath))
#all_rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", wkpath))

############################################
# Control Samples Collection for Annotations
############################################

controlSampleIDs <- anno[ (mainGroup == "A_level_control" | mainGroup == "B_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs <- controlSampleIDs[controlSampleIDs %in% colnames(all_rsemtpm)]
saveRDS(controlSampleIDs, sprintf("%s/aggregated_results/controlSampleIDs.rds", wkpath))

controlSampleIDs2 <- anno[ (mainGroup == "A_level_control") & th1 > 6000 ]$WTA.plate
controlSampleIDs2 <- controlSampleIDs2[controlSampleIDs2 %in% colnames(all_rsemtpm)]
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
#all_vars <- readRDS(sprintf("%s/aggregated_results/all_vars.rds", wkpath))

#reduces mym sample names from filename to anno list name. Currently checking all anno lists
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[x]), ignore.case = T, x = colnames(all_vars$A))
  data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
colnames(all_vars$A)[tmp[idx>0]$idx] <- tmp[idx>0]$id
colnames(all_vars$B)[tmp[idx>0]$idx] <- tmp[idx>0]$id


#phase SNPs:
# there are duplicate names from double renaming. possibly phase pre merge. 
source( sprintf('%s/scripts/phaseMyAllelesBasedOnHaplotypes.R', scriptsdir) )
mym <- phaseMyAllelesBasedOnHaplotypes(m = all_vars)
saveRDS(object = mym, file = sprintf("%s/aggregated_results/haplotypePhasedADs.rds", wkpath))

#mym <- readRDS(file = sprintf("%s/aggregated_results/haplotypePhasedADs.rds", wkpath))

#Add features to snps
ASE <- AddFeaturesandPhasing(mym, anno = anno, centromeres = centromeres)
saveRDS(object = ASE, file = sprintf("%s/aggregated_results/ASE.rds", wkpath))
#write.csv(x = ASE, file = sprintf("%s/aggregated_results/ASE.csv", wkpath))
#ASE <- readRDS(file = sprintf("%s/aggregated_results/ASE.rds", wkpath))


#Allele level data:
alleles.all <- prepareAlleleData(m = ASE, useFraction = F, controlSampleIDs = controlSampleIDs)
saveRDS(alleles.all, file = sprintf("%s/aggregated_results/alleles.all.rds", wkpath))

alleles.all.frac <- prepareAlleleData(m = ASE, useFraction = T, controlSampleIDs = controlSampleIDs)
saveRDS(alleles.all, file = sprintf("%s/aggregated_results/alleles.all.frac.rds", wkpath))


#######################
# Organize Data for ML:
#######################

#Standard set of annotations for all genes
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))

#list of preset parameters for the experiment
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))
#configs$minDetectionLevel <- 50

#To me it looks like etai's df is standardized as well
# Use the first for visuals and second for ML. Curretly ran second version.
#generate_adt_adt0_adt.na_and_nonzeros_data(dirpath = wkpath, th = 4, scriptsdir = scriptsdir, normBySd = T) 
generate_adt_adt0_adt.na_and_nonzeros_data(dirpath = wkpath, th = 4, scriptsdir = scriptsdir, normBySd = F) 

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

#Generate the arms object made of SSM, MA, & allelic CN predictions
dirpath <- wkpath
datapath <- datadir
CNdir <- sprintf("%s/CN_data", dirpath)
source(sprintf("%s/scripts/run_all_CN_predictions_for_sample_nikos.R", scriptsdir))

#ready to run "Run_ML.R" script
print("Done")


### Experimentation ###
adt <- data.table(readRDS(file = sprintf("%s/aggregated_results/adt.rds", dirpath)))
adt <- adt.default <- adt[order(seqnames, start, end)]
adt <- cbind( adt[,c(1:4)], setcolorder(adt[,-c(1:4)], order(colnames(adt[,-c(1:4)]))) )

hist(adt[,c(25)])
hist(adt_old_quantile_nocent[,c("170119_A3")])
hist(adt_old_cyclicloss[,c("170119_A3")])
hist(adt_old[,c("170119_A3")])
hist(log2(rsemtpm[,c("170119_A3")]+1), xlim = c(0,20)) 

controlSampleIDs