###################################
# Input Desired script Directories:
###################################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025f_Lane1", 1)
scriptsdir <- args[1]
wkdir <- args[2]
datadir <- args[3]
experiment <- args[4] 
threads <- args[5]
mk_wkdir <- sprintf("mkdir -p  %s", wkdir)
system(mk_wkdir)

###########
# Basic QC:
###########

source( sprintf("%s/scripts/RNAseqUtils.R", scriptsdir) )
require(data.table)
#require(snow)

anno <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir)))

rnaseqOutDir <- sprintf("%s/%s/STAR", wkdir, experiment)
expr <- collectAllSTARCountsFromFileDir(workdir = rnaseqOutDir)

qc <- data.table(expr$qc)
qc_names <- rownames(expr$qc)
counts <- data.table(expr$counts)

tmp <- rbindlist(lapply(1:length(anno$Fastq_files), function(x)  {
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[x]), ignore.case = T, x = colnames(qc))
  res <- ifelse(length(res) > 0, res, grep(pattern = sprintf("^%s", anno$Fastq_files[x]), ignore.case = T, x = colnames(qc)))
  data.table(id = unlist(anno[,c("WTA.plate")][x]),  #data.table(id = anno$WTA.plate[x], 
             idx = ifelse(length(res) > 0, res, -1)) }))

setnames(qc, tmp[idx>0]$idx, as.character(tmp[idx>0]$id))
setnames(counts, tmp[idx>0]$idx, as.character(tmp[idx>0]$id))

qc2 <- round(sweep(qc, 2, colSums(qc), FUN = "/")*100)
qcbygene <- rbind(colSums(counts>0), colSums(counts>4), colSums(counts>9))
#rownames(qcbygene) <- c("th1", "th5", "th10")

qc3 <- t(rbind(qc2, qcbygene, totalReads=colSums(qc)))
qcs <- cbind(id = rownames(qc3), qc3)
colnames(qcs) <- c("id", qc_names, "th1", "th5", "th10", "totalReads")


mk_analysis_dir <- sprintf("mkdir -p  %s/%s/QC", wkdir, experiment)
system(mk_analysis_dir)

write.csv(x = qc3, file = sprintf("%s/%s/QC/%s_QC.csv", wkdir, experiment, experiment)) #or should i write out qcs

saveRDS(counts, sprintf("%s/%s/QC/%s_counts.rds", wkdir, experiment, experiment))

#hist(qcs$geneCounts, main="", xlab="Gene Counts", cex.lab=1.5, cex.axis=1.4, col="darkgreen", xlim=c(0,100), breaks=34)


#########################################################################
## Collection of total expression data (controls and samples to predict):
#########################################################################
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

#############################################################################################################################
#Collection of SNP data (controls and samples to predict) and phasing them based on genotype data and expression (biallelic):
#############################################################################################################################

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

mk_var_dir <- sprintf("mkdir -p %s/%s/variant_results/",  wkdir,experiment)
system(mk_var_dir)
saveRDS(object = mym, file = sprintf("%s/%s/variant_results/%s_raw_variants.rds", wkdir,experiment, experiment))


