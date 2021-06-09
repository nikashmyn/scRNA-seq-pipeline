# Set path prefixs
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]

### Import Data ###

#Gene annotations
adt <- readRDS(sprintf("%s/aggregated_results/adt.rds", dirpath))
ganno <- data.table(readRDS(file = sprintf("%s/aggregated_results/ganno.rds", dirpath))) #called arms in visualization script

#P-value matrices
pval_matrix_loss_byarm <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_loss_byarm.rds", dirpath))
pval_matrix_normal_byarm <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_normal_byarm.rds", dirpath))
pval_matrix_gain_byarm <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_gain_byarm.rds", dirpath))
pval_matrix_control_byarm <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_control_byarm.rds", dirpath))

#Arm level CN Predictions from OLR
NUM_OF_FEATURES <-  50
ourpreds <- readRDS(file = sprintf("%s/ML_data/preds.probs_olr_model.WIN%d.rds", dirpath, NUM_OF_FEATURES))

#Add arm level annotations to ourpreds
preds <- data.table(ourpreds[["cns"]])
ganno2 <- ganno[,c(1,6:9)]
preds2 <- setkey(preds, seqnames, start, end, id)
ganno3 <- setkey(ganno2, seqnames, start, end, id)
preds3 <- merge(ganno3, preds2)
OLR_preds <- cbind(preds3[,c("id", "seqnames", "arm", "start", "end")], preds3[,-c(1:5)])

#total expression from OLR binned by arm
OLR_preds_byarm <- OLR_preds[,-c(1,2,4,5)] %>% 
                     group_by(arm) %>%
                       summarise_all(mean, na.rm = TRUE)

#reorder OLR preds columns like pval matrix
OLR_preds_byarm <- OLR_preds_byarm[,colnames(pval_matrix_control_byarm)]

#Thresholds for classification of CN at an arm level
a <- .05

#Classify 
CN_matrix <- c()
for (i in 1:nrow(pval_matrix_normal_byarm)) {
  row <- c()
  for (j in 2:ncol(pval_matrix_normal_byarm)) {
    switch <- T
    if (OLR_preds_byarm[i,j] >= 2.66) {
      if (pval_matrix_normal_byarm[i,..j] <= a) {
        if (pval_matrix_gain_byarm[i,..j] >= a) {
          row <- append(row, 3)
          switch <- F}}}
    if (OLR_preds_byarm[i,j] <= 1.33) {
      if (pval_matrix_normal_byarm[i,..j] <= a) {
        if (pval_matrix_loss_byarm[i,..j] >= a) {
          row <- append(row, 1)
          switch <- F}}}
    if (OLR_preds_byarm[i,j] <= 2.33) {
      if (OLR_preds_byarm[i,j] >= 1.66) {
        if (pval_matrix_normal_byarm[i,..j] >= a) {
          if (pval_matrix_gain_byarm[i,..j] <= a) {
            if (pval_matrix_loss_byarm[i,..j] <= a) {
              row <- append(row, 2)
              switch <- F}}}}}
    if (switch) {
      row <- append(row, NA)}}
  CN_matrix <- rbind(CN_matrix, row)}
CN_matrix <- cbind(OLR_preds_byarm[,c(1)], CN_matrix)
colnames(CN_matrix) <- colnames(OLR_preds_byarm)
rownames(CN_matrix) <- CN_matrix$arm
#CN_matrix <- data.frame(CN_matrix)

#get stats on CN_matrix 
#(a = .05, CN cutoffs at 1.33 1.66 2.33 and 2.66, gave frac of .57)
#(a = .05, CN cutoffs at 1.5 and 2.5, gave frac of .61)
#(a = .1, CN cutoffs at 1.5 and 2.5, gave frac of .68)
frac_confident <- length(CN_matrix[!is.na(CN_matrix)]) / (nrow(CN_matrix) * ncol(CN_matrix))

### family tree dendogram ###

#Sample families to exclude
exclude <- c("113015_T12", "061516_T12", "070516_LookRNAseq", "072516_LookRNAseq", 
             "071816_LookRNAseq", "071816_LookRNAseq", "161005_LookRNAseq", 
             "061316_LookRNAseq", "161005_LookRNAseq_control", "Jinyu", "170726_MCF10A_control",
             "040416_ATCC", "161017_LookRNAseq", "062716_LookRNAseq", "080816LookRNAseq", "080816_LookRNAseq",
             "071416_LookRNAseq", "161003_LookRNAseq")

#read in annotations
anno <- data.table(readRDS(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))
samples_to_use <- anno[!LookRNAseq.Exp %in% exclude]$WTA.plate
columns <- colnames(adt)[-c(1:4)]
samples_to_use <- c(intersect(columns, samples_to_use))

setkey(anno, WTA.plate)
families <- sort(table(anno[samples_to_use][Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
myfamily = names(families)[20]
myids <- anno[Pairs %in% myfamily]$WTA.plate
CN_matrix_byfamily <- CN_matrix[,myids]

#table of classifications
levels=c(1, 2, 3) #all unique values in df
CN_counts <- sapply(levels, function(x) rowSums(CN_matrix_byfamily==x, na.rm = T)) #count occurrences of x in each row
colnames(CN_counts) <- levels 

#Here we need to get the family information before we can get anymore accuracy with the table.



# Create a first line
plot(CN_matrix[,myids[1]], type = "l", frame = FALSE, 
     col = "red", xlab = "x", ylab = "y", ylim = c(.9, 3.1))
for (i in 2:length(myids)) {
  # Add a second line
  lines(CN_matrix[,myids[i]], col = i, type = "l")
}

histogram(CN_matrix[,4])
