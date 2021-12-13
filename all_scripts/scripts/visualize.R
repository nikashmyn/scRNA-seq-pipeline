#Visualization script CZ style 11_23_21

########################
### Passed Arguments ###
########################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
print(args)
scriptsdir <- args[1]
wkpath <- dirpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]

####################
### Read-in data ###
####################

#Read in annotation list with sister and cousin information
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
controlSampleIDs2 <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlIDs <- readRDS(file = sprintf("%s/aggregated_results/reduced_controlIDs.rds", dirpath))
arms <- ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))) #called arms in visualization script
centromeres <- readRDS( sprintf("%s/centromeres.rds", datadir) )

#raw tpm object
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))

#make new TPM object that is just raw TPM in a table with annotations
rsemtpm2 <- data.table(cbind(id = rownames(rsemtpm), rsemtpm))
setkey(rsemtpm2, "id"); setkey(geneRanges, "id")
rsemtpm3 <- merge(geneRanges, rsemtpm2)

#Read in QC object and establish high QC criteria
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))
high_qc_ids <- as.character(all_QC[which(all_QC$th5 >= quantile(all_QC$th5, c(.10))),]$id)

#Variant matrix in the coding and UTR regions
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))

#Parameters
cell <- "210111_4A"
chr <- "chr1"
bin_size <- 50

#Get TPM information for specific chr
rsemtpm_specific <- cbind(cbind(rsemtpm3[,c(1:6)], rsemtpm3[,..cell]), control_TPM = rowMeans(rsemtpm3[,..controlSampleIDs]))
rsemtpm_specific <- rsemtpm_specific[which(rsemtpm_specific$seqnames %in% chr),]

#exclude based on average ctrl tpm
rsemtpm_specific <- rsemtpm_specific[which(rsemtpm_specific$control_TPM > 5)]

#add binning to tpm object
num_of_bins <- ceiling(nrow(rsemtpm_specific)/bin_size)
rsemtpm_specific <- cbind(rsemtpm_specific, bin = rep(c(1:num_of_bins), each=bin_size)[1:nrow(rsemtpm_specific)])

#take mean by bin
rsemtpm_specific_binned <- data.table(rsemtpm_specific[,-c(1:6)] %>%
  group_by(bin) %>%
  summarise_all(mean, na.rm = TRUE))

#calculate the ratio of cell to control panel
rsemtpm_specific_binned$ratio <- (rsemtpm_specific_binned[,..cell] / rsemtpm_specific_binned$control_TPM)

#plot the tpm data across the chr
plot()
