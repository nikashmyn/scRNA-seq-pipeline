########################
### Passed Arguments ###
########################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data")
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]
mk_vis_dir <- sprintf("mkdir -p %s/visual_results", dirpath)
system(mk_vis_dir)

#####################
### Load Packages ###
#####################

require(data.table)
require(gtools)
require(readxl)
require(ggplot2)
require(matrixStats)
require(reshape2)
require(gridExtra)
require(dplyr)
require(ggthemes)

source(sprintf("%s/scripts/CNML.R", scriptsdir))
source(sprintf("%s/plots/ClusterMap_byfamily.R", scriptsdir))
source(sprintf("%s/plots/class_prob_to_cn_level.R", scriptsdir))
source(sprintf("%s/plots/plot_cn_prediction_probabilities.R", scriptsdir))
source(sprintf("%s/plots/calc_and_plot_intermediate_pvalue.R", scriptsdir))
source(sprintf("%s/plots/plot_raw_data_and_prediction_boxplots4.R", scriptsdir))
source(sprintf("%s/plots/plot_barplots_of_AllelicAndExpBiasPerSamples2.R", scriptsdir))
source(sprintf('%s/plots/plot_cn_and_intermediate_prediction_latest_model_intermediate_bars.R', scriptsdir))
source(sprintf('%s/plots/plot_pvals_and_tpm_distributions.R', scriptsdir))
source(sprintf('%s/plots/plot_pvals_and_tpm_distributions_bychr.R', scriptsdir))
source(sprintf("%s/plots/TPM_ratio_func_of_gene_expr.R", scriptsdir))

####################
### loading data ###
####################

#############################Parameters##################################

#Sample families to exclude
exclude <- c("113015_T12", "061516_T12", "070516_LookRNAseq", "072516_LookRNAseq", 
             "071816_LookRNAseq", "071816_LookRNAseq", "161005_LookRNAseq", 
             "061316_LookRNAseq", "161005_LookRNAseq_control", "Jinyu", "170726_MCF10A_control",
             "040416_ATCC", "161017_LookRNAseq", "062716_LookRNAseq", "080816LookRNAseq", "080816_LookRNAseq",
             "071416_LookRNAseq", "161003_LookRNAseq")

#Window size used in models
NUM_OF_FEATURES <- 100

##########################Data for Raw Plots###########################
message("Loading data...")

#load config file
configs <- readRDS(sprintf("%s/param_config_list.rds", datadir))

#log-space centered tpm object
adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath))
adt.na <- readRDS(sprintf("%s/aggregated_results/adt.na.rds", dirpath))
adt <- adt.default <- adt[order(seqnames, start, end)]
adt <- cbind( adt[,c(1:4)], setcolorder(adt[,-c(1:4)], order(colnames(adt[,-c(1:4)]))) )

#Create a dummy object to feed into ML preds. 
adt_fake <- adt[,-c(1:4)]
adt_fake[,] <- 0
adt_fake <- cbind(adt[,c(1:4)], adt_fake)

#raw tpm object
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))

#Read in QC object and establish high QC criteria
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))
threshold <- quantile(all_QC$th5, c(.10))
high_qc_ids <- as.character(all_QC[which(all_QC$th5 >= 4000),]$id)

#Variant matrix in the coding and UTR regions
coding <- readRDS(sprintf("%s/aggregated_results/ASE.coding.rds", dirpath))

#Low var counts exclusion. Note: both distributions of var counts are fairly similar.
#high_varqc_idsA <- names(which(colSums(coding$cnts.A[,-c(1:2)]) > quantile(colSums(coding$cnts.A[,-c(1:2)]), c(.10))))
#high_varqc_idsB <- names(which(colSums(coding$cnts.B[,-c(1:2)]) > quantile(colSums(coding$cnts.B[,-c(1:2)]), c(.10))))
#high_varqc_ids <- intersect(high_varqc_idsA, high_varqc_idsB)
#high_qc_ids <- intersect(high_qc_ids, high_varqc_ids)

#load cell annotations
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
columns <- colnames(adt)[-c(1:4)]
samples_to_use <- c(intersect(columns, anno$WTA.plate))
high_qc_ids <- intersect(high_qc_ids, samples_to_use)

#columns annotations used in file naming convention for visuals
col_anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
dim(col_anno)
excluded_by_qc <- col_anno$WTA.plate[-which(col_anno$WTA.plate %in% high_qc_ids)]
col_anno <- col_anno[ WTA.plate %in% high_qc_ids]
dim(col_anno)
col_anno[Pairs == "NA", Pairs := NA]

#fraction of nonzero tpm (binned) object
#nonzeros.zs <- readRDS(sprintf("%s/aggregated_results/nonzeros.zs.bin50.rds", dirpath))

#Gene and cell annotaions
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
controlSampleIDs2 <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs2.rds", dirpath))
controlIDs <- readRDS(file = sprintf("%s/aggregated_results/reduced_controlIDs.rds", dirpath))
arms <- ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))) #called arms in visualization script
centromeres <- readRDS( sprintf("%s/centromeres.rds", datadir) )

#Exclude this noisy control sample
#controlSampleIDs <- controlSampleIDs[-which(controlSampleIDs %in% "170425_A3")]
#controlSampleIDs2 <- controlSampleIDs2[-which(controlSampleIDs2 %in% "170425_A3")]

##########################Data for OLR Plots###########################
#message("Loading more data...")
#
##Main ML Predictions object
#ourpreds <- readRDS(file = sprintf("%s/ML_data/preds.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))
#predicted_samples <- names(ourpreds$preds)
#samples_to_use <- c(intersect(predicted_samples, samples_to_use))
#
##Auxiliary ML Data Objects
#interstat.chr <- data.table(readRDS(file = sprintf("%s/ML_data/interstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
#interstat.arm <- data.table(readRDS(file = sprintf("%s/ML_data/interstat.armlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))
#fracstat <- data.table(readRDS(file = sprintf("%s/ML_data/fracstat.chrlevel.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES)))

##########################Data for Barcharts###########################
message("Loading more data...")

allele_frac.bychr <- readRDS(file=sprintf("%s/aggregated_results/allele_frac.bychr.rds", dirpath))
chr.Af2 <- copy(allele_frac.bychr$A)
setnames(chr.Af2, old = "seqnames", new = "bin_id")
chr.Af2$bin_id <- as.character(paste0("chr", chr.Af2$bin_id))
chr.Af2$bin_id[chr.Af2$bin_id == "chr23"] <- "chrX" #Change chr23 to chrX
chr.Bf2 <- copy(allele_frac.bychr$B)
setnames(chr.Bf2, old = "seqnames", new = "bin_id")
chr.Bf2$bin_id <- as.character(paste0("chr", chr.Bf2$bin_id))
chr.Bf2$bin_id[chr.Bf2$bin_id == "chr23"] <- "chrX" #Change chr23 to chrX

#MOVE TO VISUAL SCRIPT
allele_frac.bybin <- readRDS(file=sprintf("%s/aggregated_results/allele_frac.bybin.rds", dirpath))

#########################Data for Pval Plot###########################

#calculate pvals
#source(sprintf("%s/scripts/calculating_pvals_byarm_nikos.R", scriptsdir))
source(sprintf("%s/scripts/pval_grouped_samples.R", scriptsdir))
source(sprintf("%s/scripts/calculating_control_pvals_bychr_v3.R", scriptsdir))
source(sprintf("%s/scripts/calculating_pvals_byarm_v2.R", scriptsdir))

#Add arm level annotations to ourpreds
#preds <- data.table(ourpreds[["cns"]])
#ganno2 <- ganno[,c(1,6:9)]
#preds2 <- setkey(preds, seqnames, start, end, id)
#ganno3 <- setkey(ganno2, seqnames, start, end, id)
#preds3 <- merge(ganno3, preds2)
#OLR_preds_warm <- cbind(preds3[,c("arm")], preds3[,-c(1:5)])
#
##total expression from OLR binned by arm
#OLR_preds_byarm <- OLR_preds_warm %>% 
#  group_by(arm) %>%
#  summarise_all(mean, na.rm = TRUE)

#reorder OLR preds columns like pval matrix
#OLR_preds_byarm <- OLR_preds_byarm[,colnames(pval_matrix_control_byarm)]


########################Data for Cluster Plot##########################

#TPM object normalized with chr bins
adt_bychr <- readRDS(file=sprintf("%s/aggregated_results/adt.bychr.rds", dirpath))
#absolute difference between allele cnts per chr
abs_allele_diff_mat.bychr <- readRDS(file=sprintf("%s/aggregated_results/abs_allele_diff_mat.bychr.rds", dirpath))
#total SNP counts over both alleles taken after separate allele normalization. ~same as var_mat.bychr. 
var_mat_sep.bychr <- readRDS(file=sprintf("%s/aggregated_results/var_mat_sep.bychr.rds", dirpath))

adt_bychr2 <- cbind( bin_id = paste0("chr", unlist(adt_bychr[,1])), adt_bychr[,-1]*2)
adt_bychr2$bin_id[adt_bychr2$bin_id=="chr23"] <- "chrX"

#add in quadrature and scale by sqrt(2) so that the point (1,1) is 1 in aggregated dim
matrix_adt.bychr <- as.matrix(adt.bychr[,-c(1)]); matrix_var.bychr <- as.matrix(var_mat_sep.bychr[,-c(1)]);
matrix_CN_Ratio <- sqrt( ( (matrix_adt.bychr^2) + (matrix_var.bychr^2) ) ) / sqrt(2)
CN_ratio.bychr <- data.table(cbind(adt.bychr[,c(1)], matrix_CN_Ratio))

#Read in Pvals
non_diploid_chrs_w_vals <- readRDS(file = sprintf("%s/aggregated_results/non_diploid_chromosomes.rds", dirpath))

message("Loading completed!")


##################################
### Read in possible MN events ###
##################################

possible_events <- read.csv(file = sprintf("%s/aggregated_results/possible_MN_events.csv", dirpath))
no_events <- read.csv(file = sprintf("%s/aggregated_results/cells_noMN_events.csv", dirpath))
all_family_events <- data.table(read.csv(file = sprintf("%s/aggregated_results/all_family_events.csv", dirpath)))

############################
### adt no log all genes ###
############################

dt <- data.table(rsemtpm, keep.rownames = "id")
if(length(names(which(table(colnames(dt))>1))) > 0)
  dt <- dt[, -which(colnames(dt) %in% names(which(table(colnames(dt))>1))), with=F]
setkey(geneRanges, id)
setkey(dt, id)
dt <- merge(geneRanges, dt)
require(gtools)
dt$seqnames <- as.character(dt$seqnames)
dt$seqnames <- factor(dt$seqnames, mixedsort(unique(dt$seqnames)))
dt <- dt[order(seqnames, start, end)]
adt_nolog <- dt[seqnames != "chrM" & seqnames != "chrY"]

write.csv(adt_nolog, file = sprintf("%s/aggregated_results/raw_tpm.csv", dirpath), row.names=FALSE)

#####################################
### Run Visual Report in Parallel ###
#####################################

require(doParallel)
registerDoParallel(20)

#Pull information for a given sample from ML prediction object
#get_flat_format_tbl <- function(mysampid) {
#  probs <- data.table(sample_id = mysampid, cbind(ourpreds$rowdata,  ourpreds$preds[[mysampid]]))
#  setkeyv(probs, c("seqnames", "start", "end"))
#  setkeyv(arms, c("seqnames", "start", "end"))
#  tmp <- foverlaps(probs, arms)
#  tmp2 <- tmp[, c("sample_id", "seqnames", "arm", "i.start", "i.end", "id", "1", "2", "3"), with=F]
#  setnames(x = tmp2, old = c("i.start", "i.end", "1", "2", "3"), new = c("start", "end", "c1", "c2", "c3"))
#  return(tmp2)
#}


#Variables
destDir <- sprintf("%s/visual_results", dirpath)
system(sprintf("mkdir -p %s", destDir))
#chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr11", "chr12", "chr16", "chr17", "chr19")
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")
chrs <- chrs[which(!chrs %in% configs$chr_to_excl)]
#arms_to_exclude <- c("10q", "13q", 18", "20", "21", "22", "X")
# 14, 15 have had no p arm values and the rest have intentionally excluded arms 

#Get the samples that can be visualized
samples_to_visualize <- intersect(col_anno$WTA.plate, samples_to_use)
exclude_chrs <- c()

#make raw plot pdfs for every sample
foreach(i = c(1:length(samples_to_visualize))) %dopar% {
  myid <- samples_to_visualize[i]
  #probs <- get_flat_format_tbl(myid)
  myfamily_file <- anno[which(anno$WTA.plate == myid),]$Family_IDs
  #if(myfamily == "") {myfamily <- "no_family"}
  #myfamily_file1 = str_remove(myfamily, fixed("(")); myfamily_file = str_remove(myfamily_file1, fixed(")"))
  myid_file <- anno[which(anno$WTA.plate == myid),]$Cell_IDs
  message(myfamily_file)
  #Plot raw plots for relevant chromosomes
  for(chr in chrs) {
    #Chr information
    message(chr)
    #Plot raw data visuals function
    plot_raw_data_and_prediction_boxplots4(myfamily_file = myfamily_file, myid_file = myid_file, myid = myid, chr = chr, adt = adt, destDir = destDir,
                                           coding = coding, th=configs$minNumOfSamplesToDetect, 
                                           preds = adt_fake, controlSampleIDs = controlSampleIDs)
    CumulativeTPM_and_ratioTPM_plots(adt = adt_nolog, geneRanges = geneRanges, controlSampleIDs = controlSampleIDs, 
                                     myfamily_file = myfamily_file, myid_file = myid_file, myid = myid, chr = chr, destDir = destDir, binsize = 10)
  }
}


#Make barcharts by family
setkey(anno, WTA.plate)
families <- sort(table(anno[samples_to_visualize][Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }
exclude_chrs <- c()

make_path <- sprintf("mkdir -p %s/byfamily/", destDir)
system(make_path)

foreach(i = c(1:length(names(families)))) %dopar% {
#for(i in 1:2) {
  myfamily = names(families)[i]
  myids <- anno[Pairs %in% myfamily]$WTA.plate
  myids_name <- anno[Pairs %in% myfamily]$Cell_IDs
  #myfamily_file1 = str_remove(myfamily, fixed("(")); myfamily_file = str_remove(myfamily_file1, fixed(")"))
  fam_id_name <- unique(anno$Family_IDs[which(anno$Pairs == myfamily)])
  message(myfamily) 
  #plot Barcharts by family
  plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = adt_bychr2, myfamily_file = fam_id_name, destDir = sprintf("%s/byfamily", destDir),
                                               wt = rep(1:nrow(adt_bychr2)), fracs = list(Af = chr.Af2, Bf = chr.Bf2), nogeno = c(),
                                               ids = myids, name_ids = myids_name, chr = "", plotAxisText = F)
}

#Reduce anno list to relevant gen2 samples
rows1 <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
rows2 <- intersect(which(anno$key_pairs == "Gen1"), which(anno$exclusion == ""))
rows <- union(rows1, rows2)
anno_gen_1_2 <- anno[rows,]

#extract family information
setkey(anno_gen_1_2, WTA.plate)
families <- sort(table(anno_gen_1_2[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }

make_path <- sprintf("mkdir -p %s/byfamily/", destDir)
system(make_path)

foreach(i = c(1:length(names(families)))) %dopar% {
  #for(i in 1:2) {
  myfamily = names(families)[i]
  myfamily_file1 = str_remove(myfamily, fixed("(")); myfamily_file = str_remove(myfamily_file1, fixed(")"))
  message(myfamily_file)
  myids <- anno[Pairs %in% myfamily]$WTA.plate
  pdf(file = sprintf("%s/byfamily/%s.family.pdf", destDir, myfamily_file), width = 18, height = 12)
  #get events for this family
  possible_events_myfam <- all_family_events[which(all_family_events$family == myfamily),]
  #no_events_myfam <- no_events[which(no_events$family == myfamily),]
  #if(nrow(possible_events_myfam) == 0) {
  #  if(nrow(no_events_myfam)>0) {
  #    possible_events_myfam <- no_events_myfam[which(no_events_myfam$chr == 1),]}}
  #plot table of possible events
  if (nrow(possible_events_myfam) >= 1) {grid.table(possible_events_myfam)}
  #plot Barcharts by family
  q <- plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = adt_bychr2, myfamily_file = myfamily_file, destDir = sprintf("%s/byfamily", destDir),
                                                    wt = rep(1:nrow(adt_bychr2)), fracs = list(Af = chr.Af2, Bf = chr.Bf2), nogeno = c(),
                                                    ids = myids, chr = "", plotAxisText = F, return_plot = T)
  #Plot cluster map output by family
  #ClusterMap2D_byfamily(adt.bychr, var_mat_sep.bychr, myids, myfamily, exclude_chrs)
  #plot information for possible events
  if (nrow(possible_events_myfam) >= 1 ) {
    #Plot raw plots for relevant chromosomes
    for(chr in unique(possible_events_myfam$chr)) {
      #Chr information
      chr <- paste0("chr", chr)
      if(chr == "chr23") {chr <- "chrX"}
      message(chr)
      files <- data.table()
      files2 <- data.table()
      for(myid in myids[myids %in% samples_to_visualize]) {
        #Plot raw data visuals function
        files <- rbind(files, plot_raw_data_and_prediction_boxplots4(myfamily_file = myfamily_file, myid = myid, chr = chr, adt = adt, destDir = destDir,
                                                                     coding = coding, th=configs$minNumOfSamplesToDetect,
                                                                     preds = adt_fake, controlSampleIDs = controlSampleIDs, save.rds = T))
        files2 <- rbind(files2, CumulativeTPM_and_ratioTPM_plots(adt = adt_nolog, geneRanges = geneRanges, controlSampleIDs = controlSampleIDs, save.rds = T,
                                                                 myfamily_file = myfamily_file, myid = myid, chr = chr, destDir = destDir, binsize = 10))
      }
      #TPMratio by bin and var by bin raw plots
      TPM_vis <- lapply(files$TPM_file , function(x){readRDS(x)} )
      var_vis <- lapply(files$var_file , function(x){readRDS(x)} )
      interweave <- c()
      for(i in 1:length(TPM_vis)){interweave <- append(interweave, TPM_vis[i]); interweave <- append(interweave, var_vis[i])}
      do.call("grid.arrange", c(interweave, nrow = length(TPM_vis), ncol = 2)) 
      #CumulativeTPM and ratio_to_expr plots
      expr_vis <- lapply(files2$expr_file , function(x){readRDS(x)} )
      cum_vis <- lapply(files2$cum_file , function(x){readRDS(x)} )
      interweave <- c()
      for(i in 1:length(expr_vis)){interweave <- append(interweave, expr_vis[i]); interweave <- append(interweave, cum_vis[i])}
      do.call("grid.arrange", c(interweave, nrow = length(expr_vis), ncol = 2)) 
    }
  }
  dev.off()
}

#Sanity check that the code ran to completion
print("Done with making visuals")


########################


df <- cbind(rep("Intact", 8), rep("No Defect", 8))
df <- rbind(df, cbind(rep("Early NE Defect", 10), rep("No Defect", 10)))
df <- rbind(df, cbind(rep("S/G2", 10), rep("No Defect", 10)))
df <- rbind(df, cbind(rep("Early NE Defect", 4), rep("Defect", 4)))
df <- rbind(df, cbind(rep("S/G2", 8), rep("Defect", 8)))
colnames(df) <- c("rupt_time", "state")
df <- df[,c("state", "rupt_time")]
df <- data.table(df)

df2 <- cbind(rep("Intact", 4), rep("No Defect", 4))
df2 <- rbind(df2, cbind(rep("Early NE Defect", 1), rep("No Defect", 1)))
df2 <- rbind(df2, cbind(rep("Intact", 8), rep("Defect", 8)))
df2 <- rbind(df2, cbind(rep("Early NE Defect", 8), rep("Defect", 8)))
df2 <- rbind(df2, cbind(rep("S/G2", 6), rep("Defect", 6)))
colnames(df2) <- c("rupt_time", "state")
df2 <- df2[,c("state", "rupt_time")]
df2 <- data.table(df2)



p <- ggplot(df, aes(x = rupt_time, fill = state )) + 
  geom_bar() + 
  scale_fill_manual(values = c("black", "grey")) +
  #scale_y_discrete(breaks = c(0,5,10,15,20)) +
  scale_y_continuous(breaks = c(0,4,8,12,16,20), labels = c("0","4","8","12","16","20"), limits = c(0,20)) + 
  scale_x_discrete(limits = c("Intact", "Early NE Defect", "S/G2")) +
  coord_cartesian(clip="off") +
  geom_rangeframe(x=0, y=seq(0, 20, along.with = df$state)) + # 
  theme_tufte() +
  theme(aspect.ratio = 2, text=element_text(size=18), axis.ticks.x = element_blank(),
        axis.text = element_text(color="black"), axis.title = element_blank()) 


p1 <- ggplot(df2, aes(x = rupt_time, fill = state )) + 
  geom_bar() + 
  scale_fill_manual(values = c("black", "grey")) +
  #scale_y_discrete(breaks = c(0,5,10,15,20)) +
  scale_y_continuous(breaks = c(0,4,8,12), labels = c("0","4","8","12"), limits = c(0,12)) + 
  scale_x_discrete(limits = c("Intact", "Early NE Defect", "S/G2")) +
  coord_cartesian(clip="off") +
  geom_rangeframe(x=0, y=seq(0, 12, along.with = df2$state)) + # 
  theme_tufte() +
  theme(aspect.ratio = 2, text=element_text(size=18), axis.ticks.x = element_blank(),
        axis.text = element_text(color="black"), axis.title = element_blank()) 

ggsave(plot = p, filename = sprintf("Gen2_barchart.pdf"), path = destDir, device = "pdf", width = 6, height = 12, dpi = 300, units = "in")
ggsave(plot = p1, filename = sprintf("Gen1_barchart.pdf"), path = destDir, device = "pdf", width = 6, height = 12, dpi = 300, units = "in")


