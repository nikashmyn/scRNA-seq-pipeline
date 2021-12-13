##########################
### Script Explanation ###
##########################

#This is the first analysis scipt that is run as a part of our scRNA-seq snakemake pipeline 
#This script aggregates all the alignment, expression, variant information for each seperate experiment. 
#This script feeds into "Data_Aggregation.R" which further aggregates this data. 

#-------------------------------------------------------
#This script was written by Nikos Mynhier and Etai Jacob.

########################################
### Input Desired script Directories ###
########################################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a", 1)
scriptsdir <- args[1]
wkdir <- args[2]
datadir <- args[3]
experiment <- args[4] 
threads <- args[5]
mk_wkdir <- sprintf("mkdir -p  %s", wkdir)
system(mk_wkdir)

####################################
### Collect Basic QC Information ###
####################################

source( sprintf("%s/scripts/RNAseqUtils.R", scriptsdir) )
require(data.table)

#read in cell annotations
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))

#Collect QC information post alignement from STAR output files
rnaseqOutDir <- sprintf("%s/%s/STAR", wkdir, experiment)
expr <- collectAllSTARCountsFromFileDir(workdir = rnaseqOutDir)

#Reorganize the QC information
qc <- data.table(expr$qc)
qc_names <- rownames(expr$qc)
counts <- data.table(expr$counts)

#IDs are currently named after STAR files. Change to sample ids from cell annotations. 
tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[x]), ignore.case = T, x = colnames(qc))
  res <- ifelse(length(res) > 0, res, grep(pattern = sprintf("^%s", anno$Fastq_files[x]), ignore.case = T, x = colnames(qc)))
  data.table(id = unlist(anno[,c("WTA.plate")][x]),  #data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))
setnames(qc, tmp[idx>0]$idx, as.character(tmp[idx>0]$id))
setnames(counts, tmp[idx>0]$idx, as.character(tmp[idx>0]$id))

#Divide by cols by col sums and multiply by 100 to get a percent
qc2 <- round(sweep(qc, 2, colSums(qc), FUN = "/")*100)
#Calculate th1, th5, th10. Genes covered by at least 1, 5, or 10 reads.
qcbygene <- rbind(colSums(counts>0), colSums(counts>4), colSums(counts>9))

#Restructure information and label columns
qc3 <- t(rbind(qc2, qcbygene, totalReads=colSums(qc)))
qcs <- cbind(id = rownames(qc3), qc3)
colnames(qcs) <- c("id", qc_names, "th1", "th5", "th10", "totalReads")

#make analysis directory based on input file path
mk_analysis_dir <- sprintf("mkdir -p  %s/%s/QC", wkdir, experiment)
system(mk_analysis_dir)

#Write QC information for this experiment to new file path
write.csv(x = qc3, file = sprintf("%s/%s/QC/%s_QC.csv", wkdir, experiment, experiment)) #or should i write out qcs

#Save count information as well. 
saveRDS(counts, sprintf("%s/%s/QC/%s_counts.rds", wkdir, experiment, experiment))

#############################################################################
### Collection of total expression data (controls and samples to predict) ###
#############################################################################
myrsem <- collectAllRSEMresultsFromDir(dir = sprintf("%s/%s/RSEM/output", wkdir, experiment))
mk_expr_dir <- sprintf("mkdir -p %s/%s/expr_results", wkdir, experiment)
system(mk_expr_dir)
saveRDS(myrsem, sprintf("%s/%s/expr_results/%s_rsemoutput.rds", wkdir, experiment, experiment))

myrsemtpm <- myrsem$genes$TPM

mytmpids <- do.call(rbind, lapply(rownames(myrsemtpm), function(x) str_split(string = x, pattern = "_", n = 2)[[1]]))[,1]
mytmpids[grep("_PAR_Y", x = rownames(myrsemtpm))] <- paste(mytmpids[grep("_PAR_Y", x = rownames(myrsemtpm))], "_PAR_Y", sep = "")
rownames(myrsemtpm) <- mytmpids

saveRDS(myrsemtpm, sprintf("%s/%s/expr_results/%s_rsemtpm.rds", wkdir, experiment, experiment))
write.csv(x = myrsemtpm, file = sprintf("%s/%s/expr_results/%s_rsemtpm.csv", wkdir, experiment, experiment))

###################################################################################################################################
### Collection of SNP data (controls and samples to predict) and phasing them based on genotype data and expression (biallelic) ###
###################################################################################################################################

#Parameters for Genotyping script
#--------------------------------
source( sprintf("%s/scripts/runRNAseqPipelineFor_BF9VP.utils.R", scriptsdir) )
wkdir_cur <- sprintf("%s/%s", wkdir, experiment)
numOfCluster <- as.numeric(threads) 
outpath <- sprintf("%s/Variants", wkdir_cur)
patterns <- sprintf("%s/Variants/%s_RPE_hets.GT_UG.%s.vcf", wkdir_cur, experiment, '%s')
source(sprintf("%s/scripts/prepareGenotypeDataFromRNAvcfs.R", scriptsdir))

#Load and phase data from specific lib:
#--------------------------------------
source( sprintf("%s/scripts/collectAllVCFdata.R", scriptsdir) )
alleles_AD_files <- sprintf("%s/%s/Variants/alleles.ADs.v5.chr%s.rds", wkdir, experiment, '%s')
mym <- collectAllVCFdata(patternFile = alleles_AD_files)

#Create new directory and save data
#----------------------------------
mk_var_dir <- sprintf("mkdir -p %s/%s/variant_results/",  wkdir,experiment)
system(mk_var_dir)
saveRDS(object = mym, file = sprintf("%s/%s/variant_results/%s_raw_variants.rds", wkdir,experiment, experiment))


