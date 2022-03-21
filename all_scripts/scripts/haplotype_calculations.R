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
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
print(args)
scriptsdir <- args[1]
wkpath <- dirpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]

########################################
### Read in Hets file and Phase file ###
########################################

#read in Anno list
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))

#read in phase from Greg Brunette
phase_file <- "/pellmanlab/stam_niko/test_ASE/Alleles_v2/RPE1_Haplotype_update.dat"
hets_phase <- data.table(read_table(phase_file))
hets_phase <- hets_phase[,c("chr", "pos", "Haplotype")]
setkey(hets_phase, "chr", "pos", "Haplotype")

#Generate the directories based on experiments
dir_paths <- list()
for (i in experiments) {dir_paths <- cbind(dir_paths, sprintf("%s/%s/ASE", wkpath, i))}

#list all files in each directory folder
files <- c()
for (i in 1:length(dir_paths)){files <- append(files, paste(dir_paths[i], list.files(path = as.character(dir_paths[i]), pattern = ".table"), sep = "/"))}

#List sample names
files_samples <- files
for (i in 1:length(anno$Fastq_files)){
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[i]), ignore.case = T, value=T, x = files)
  if (length(res) < 1) {res <- grep(pattern =  anno$Fastq_files[i], ignore.case = T, value=T, x = files)}
  files_samples[which(files == res[1])] <- anno$WTA.plate[i]
}
files_samples_v2 <- files_samples[which(files_samples %in% anno$WTA.plate)]
files_v2 <- files[which(files_samples %in% anno$WTA.plate)]

#overlap phase file entries with ASE tables
hets_phase_w_counts <- hets_phase; j=0
for (i in 1:length(files_v2)) {
  #read in ASEReadCounter tables
  ASE_file <- files_v2[i]
  ASE <- data.table(read.table(ASE_file, header = T))
  if (nrow(ASE) > 0) {
    setnames(ASE, old = c("contig", "position"), c("chr", "pos"))
    ASE <- ASE[,c("chr", "pos", "refCount")]
    setkey(ASE, "chr", "pos")
    ref_counts <- merge(hets_phase, ASE, all.x = T)[,c("refCount")]
    hets_phase_w_counts <- cbind(hets_phase_w_counts, ref_counts)
    colnames(hets_phase_w_counts)[i+3-j] <- files_samples_v2[i]
  } 
  else {
    j = j + 1
    sprintf("Empty file: %s", ASE_file)
  }
}

#remove completely uncovered sites
hets_phase_w_counts_2 <- hets_phase_w_counts[which(rowSums(hets_phase_w_counts[,-c(1:3)], na.rm = T)>0)]

#Keep only heterozygous sites
hets_phase_w_counts_het <- hets_phase_w_counts_2[which(!hets_phase_w_counts_2$Haplotype == 0),]
hets_phase_w_counts_hom <- hets_phase_w_counts_2[which(hets_phase_w_counts_2$Haplotype == 0),]

#Get SNP counts by haplotype
hets_phase_w_counts_A <- hets_phase_w_counts_2[which(hets_phase_w_counts_2$Haplotype == 1),]
hets_phase_w_counts_B <- hets_phase_w_counts_2[which(hets_phase_w_counts_2$Haplotype == -1),]

#overlap geneRanges annotations
ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))) #called arms in visualization script
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
setkey(geneRanges, "seqnames", "start", "end")
#add geneRanges anno to A allele and remove sites not inside of known genes
setnames(hets_phase_w_counts_A, new = c("seqnames", "start"), old = c("chr", "pos"))
hets_phase_w_counts_A$end <- hets_phase_w_counts_A$start
setkey(hets_phase_w_counts_A, "seqnames", "start", "end")
phased_het_counts_A_bySNP <- foverlaps(hets_phase_w_counts_A, geneRanges, by.x = c("seqnames", "start", "end"), by.y = c("seqnames", "start", "end"))
setnames(phased_het_counts_A_bySNP, old = c("start", "end", "i.start", "i.end"), new = c("gene_start", "gene_end", "SNP_start", "SNP_end"))
phased_het_counts_A_bySNP_2 <- cbind(phased_het_counts_A_bySNP[,c("seqnames", "id", "gene_start", "gene_end", "width", "strand", "SNP_start", "SNP_end")], phased_het_counts_A_bySNP[,-c("seqnames", "id", "gene_start", "gene_end", "width", "strand", "SNP_start", "SNP_end")])
phased_het_counts_A_bySNP_3 <- phased_het_counts_A_bySNP_2[which(!is.na(phased_het_counts_A_bySNP_2$id)),]
#add geneRanges anno to B allele and remove sites not inside of known genes
setnames(hets_phase_w_counts_B, new = c("seqnames", "start"), old = c("chr", "pos"))
hets_phase_w_counts_B$end <- hets_phase_w_counts_B$start
setkey(hets_phase_w_counts_B, "seqnames", "start", "end")
phased_het_counts_B_bySNP <- foverlaps(hets_phase_w_counts_B, geneRanges, by.x = c("seqnames", "start", "end"), by.y = c("seqnames", "start", "end"))
setnames(phased_het_counts_B_bySNP, old = c("start", "end", "i.start", "i.end"), new = c("gene_start", "gene_end", "SNP_start", "SNP_end"))
phased_het_counts_B_bySNP_2 <- cbind(phased_het_counts_B_bySNP[,c("seqnames", "id", "gene_start", "gene_end", "width", "strand", "SNP_start", "SNP_end")], phased_het_counts_B_bySNP[,-c("seqnames", "id", "gene_start", "gene_end", "width", "strand", "SNP_start", "SNP_end")])
phased_het_counts_B_bySNP_3 <- phased_het_counts_B_bySNP_2[which(!is.na(phased_het_counts_B_bySNP_2$id)),]

#Aggregate haplotype A coverage by gene
varA_mat.bygene <- phased_het_counts_A_bySNP_3[,-c(1,3:9)] %>%
  group_by(id) %>%
  summarise_all(mean, na.rm = TRUE)
varA_mat <- merge(geneRanges, varA_mat.bygene, by = "id")
varA_mat[varA_mat == "NaN"] <- 0

#Aggregate haplotype B coverage by gene
varB_mat.bygene <- phased_het_counts_B_bySNP_3[,-c(1,3:9)] %>%
  group_by(id) %>%
  summarise_all(mean, na.rm = TRUE)
varB_mat <- merge(geneRanges, varB_mat.bygene, by = "id") #merge with ganno if you want to get same genes as transcript analysis
varB_mat[varB_mat == "NaN"] <- 0

#limit to common genes
common_genes <- intersect(varA_mat$id, varB_mat$id)
varA_mat_2 <- varA_mat[which(varA_mat$id %in% common_genes)]
varB_mat_2 <- varB_mat[which(varB_mat$id %in% common_genes)]

#combine matrices
var_mat <- list(A = varA_mat_2, B = varB_mat_2)
var_mat$A <- var_mat$A[order(seqnames, start),]
var_mat$B <- var_mat$B[order(seqnames, start),]


#bin in groups of X genes
binsize <- 50
bins <- lapply(unique(var_mat$A$seqnames), function(chr) unlist(lapply(1:ceiling(nrow(var_mat$A[seqnames==chr])/binsize), rep, binsize))[1:nrow(var_mat$A[seqnames==chr])])
names(bins) <- unique(var_mat$A$seqnames)
bins2 <- rbindlist(lapply(names(bins), function(chr) data.table(seqnames = chr, bin = bins[[chr]])))
stopifnot(sum(bins2$seqnames != var_mat$A$seqnames)==0)
var_mat$A_bin <- cbind(bin = bins2$bin, var_mat$A)
var_mat$B_bin <- cbind(bin = bins2$bin, var_mat$B)

#aggregate to bin level
var_mat$A_binned <- var_mat$A_bin[,-c(2,4:7)] %>%
  group_by(seqnames, bin) %>%
  summarise_all(mean, na.rm = TRUE)

#aggregate to bin level
var_mat$B_binned <- var_mat$B_bin[,-c(2,4:7)] %>%
  group_by(seqnames, bin) %>%
  summarise_all(mean, na.rm = TRUE)

### Summary of Variant Matrix ###
#var_mat is SNP coverage averaged across gene and then averaged across bins
#SNPs were limited to be inside of a gene
#genes are limited to those with known SNPs on both haplotypes, common genes
#turned uncovered (NA) sites to 0 coverage for kept genes
saveRDS(var_mat, file=sprintf("%s/aggregated_results/ASE_matrix.rds", dirpath))
var_mat <- readRDS(file=sprintf("%s/aggregated_results/ASE_matrix.rds", dirpath))

#I need to treat the matrix how I treated the previous one. 
#look at visual code to see which one is used. 

#Normalize varA SNP ratio dataframe 
dt2go <- data.table(copy(var_mat$A_binned[, -c(1,2)]))

#Divide by row means
varA_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
varA_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), varA_rowmeans)))
dt3go <- cbind(var_mat$A_binned[, c(1,2)], (dt2go / varA_rowmeans_mat_reciprocal))
#use ratios to get cell wide bias
varA_colmeans <- colMeans(dt3go[,-c(1,2)], na.rm = T)
dtgo_biased <- data.table(1/t(replicate(nrow(dt3go), varA_colmeans)))
dtgo_biased[dtgo_biased == Inf] <- 0 
#adjust average counts by cell wide bias
dt5go_A <- cbind(var_mat$A_binned[,c(1,2)] , (dt2go * dtgo_biased))

#Normalize varB SNP ratio dataframe 
dt2go <- data.table(copy(var_mat$B_binned[, -c(1,2)]))

#Divide by row means
varB_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
varB_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), varB_rowmeans)))
dt3go <- cbind(var_mat$B_binned[, c(1,2)], (dt2go / varB_rowmeans_mat_reciprocal))
#use ratios to get cell wide bias
varB_colmeans <- colMeans(dt3go[,-c(1,2)], na.rm = T)
dtgo_biased <- data.table(1/t(replicate(nrow(dt3go), varB_colmeans)))
dtgo_biased[dtgo_biased == Inf] <- 0 
#adjust average counts by cell wide bias
dt5go_B <- cbind(var_mat$B_binned[,c(1,2)] , (dt2go * dtgo_biased))


#calculate allelic fraction by dividing each allele by the sum of the two alleles.
allele_A <- data.table(cbind(dt5go_A[,c(1,2)], (dt5go_A[,-c(1,2)])/(dt5go_A[,-c(1,2)] + dt5go_B[,-c(1,2)]) ))
allele_B <- data.table(cbind(dt5go_B[,c(1,2)], (dt5go_B[,-c(1,2)])/(dt5go_A[,-c(1,2)] + dt5go_B[,-c(1,2)]) ))
allele_frac_mat.bychr <- list(A = allele_A, B = allele_B)

#The allelic fractions are not very even for the diploid cells so limiting to common genes is distorting the balance of the SNPs. 





#Normalize varA SNP ratio dataframe 
dt2go <- data.table(copy(var_mat$A_binned[, -c(1,2)]))

#Divide by row means
varA_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
varA_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), varA_rowmeans)))
dt3go <- (dt2go / varA_rowmeans_mat_reciprocal)
#use ratios to get cell wide bias
varA_colmeans <- colMeans(dt3go, na.rm = T)
dtgo_biased <- data.table(1/t(replicate(nrow(dt3go), varA_colmeans)))
dtgo_biased[dtgo_biased == Inf] <- 0 
ratios_A <- cbind(var_mat$A_binned[,c(1,2)], (dt3go * dtgo_biased))
#adjust average counts by cell wide bias



#Normalize varB SNP ratio dataframe 
dt2go <- data.table(copy(var_mat$B_binned[, -c(1,2)]))

#Divide by row means
varB_rowmeans <- rowMeans(dt2go[,..controlSampleIDs], na.rm = T)
varB_rowmeans_mat_reciprocal <- data.table((replicate(ncol(dt2go), varB_rowmeans)))
dt3go <- (dt2go / varB_rowmeans_mat_reciprocal)
#use ratios to get cell wide bias
varB_colmeans <- colMeans(dt3go, na.rm = T)
dtgo_biased <- data.table(1/t(replicate(nrow(dt3go), varB_colmeans)))
dtgo_biased[dtgo_biased == Inf] <- 0 
ratios_B <- cbind(var_mat$B_binned[,c(1,2)], (dt3go * dtgo_biased))
#adjust average counts by cell wide bias



#calculate allelic fraction by dividing each allele by the sum of the two alleles.
var_ratios.bybin <- list(A = ratios_A, B = ratios_B)




#Plot ratios
hist(unlist(ratios_A[c(2),-c(1,2)]), breaks = 100)
abline(v = median(unlist(ratios_A[c(2),-c(1,2)])))
abline(v = mean(unlist(ratios_A[c(2),-c(1,2)])))

