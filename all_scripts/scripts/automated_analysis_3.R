#automated analysis script

##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

####################
### Read-in data ###
####################

#Read in data
ref_TPM <- readRDS(file = sprintf("%s/aggregated_results/ref_TPM.rds", dirpath))
ref_AS_TPM <- readRDS(file = sprintf("%s/aggregated_results/ref_AS_TPM.rds", dirpath))

#Read in grouped sample information
hand_samples <- readRDS(sprintf("%s/aggregated_results/grouped_control_aneuploidies.rds", dirpath))

#read in normed TPM by chr
TPM_bychr <- readRDS(file=sprintf("%s/aggregated_results/TPM.inv_var.bychr.rds", dirpath))

#new ASE
ASE_bychr <- readRDS(file=sprintf("%s/aggregated_results/ASE.inv_var.bychr.rds", dirpath))

#new AS-TPM
AS_TPM_bychr <- readRDS(file=sprintf("%s/aggregated_results/AS-TPM.inv_var.bychr.rds", dirpath))
VAR_A <- copy(AS_TPM_bychr$A)
VAR_B <- copy(AS_TPM_bychr$B)

#read in pvalues from control dist bychr
pval_matrix_control_bychr <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_control_bychr.rds", dirpath))
pval_matrix_gain_bychr <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_gain_bychr.rds", dirpath))
pval_matrix_loss_bychr <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_loss_bychr.rds", dirpath))

#read in AS_TPM pvals
AS_TPM_bychr_mono_pvals <- readRDS(file = sprintf("%s/aggregated_results/AS_TPM_pval_matrix_loss_bychr.rds", dirpath))
AS_TPM_bychr_ctrl_pvals <- readRDS(file = sprintf("%s/aggregated_results/AS_TPM_pval_matrix_control_bychr.rds", dirpath))
AS_TPM_bychr_tri_pvals <- readRDS(file = sprintf("%s/aggregated_results/AS_TPM_pval_matrix_gain_bychr.rds", dirpath))
AS_TPM_bychr_ctrl_pvals_chr_spec <- readRDS(file = sprintf("%s/aggregated_results/pval_AS_TPM_bychr_ctrl_pvals_chr_spec.rds", dirpath))
AS_TPM_bychr_mono_pvals_chr_spec <- readRDS(file = sprintf("%s/aggregated_results/pval_AS_TPM_bychr_mono_pvals_chr_spec.rds", dirpath))

#Read in annotation list with sister and cousin information
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
changeCols <- colnames(anno)[-c(13:16)]
anno[,(changeCols):= lapply(.SD, as.character), .SDcols = changeCols] 

####################
### Boxplot Data ###
####################

#store high and low cutoffs for later
low_lim <- quantile(ref_TPM$vals[which(ref_TPM$CN == 1)], .05)
high_lim <- quantile(ref_TPM$vals[which(ref_TPM$CN == 3)], .95)
low_lim_AS_TPM <- quantile(ref_AS_TPM$value[which(ref_AS_TPM$CN == 0)], .05)
high_lim_AS_TPM <- quantile(ref_AS_TPM$value[which(ref_AS_TPM$CN == 2)], .95)

#######################
### Strip Plot Data ###
#######################
#Reduce anno list to relevant gen2 samples
rows1 <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
rows2 <- intersect(which(anno$key_pairs == "Gen1"), which(anno$exclusion == ""))
rows <- union(rows1, rows2)
anno_gen_1_2 <- anno[rows,]

#extract family information
setkey(anno_gen_1_2, WTA.plate)
families <- sort(table(anno_gen_1_2[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }
exclude_chrs <- c()

names <- anno_gen_1_2$WTA.plate
ids <- anno_gen_1_2$Cell_IDs
relationship <- anno_gen_1_2$Relationship
families <- anno_gen_1_2$Pairs
family_ids <- anno_gen_1_2$Family_IDs
gen <- anno_gen_1_2$key_pairs
rupt_time <- anno_gen_1_2$MN_rupt_time_simple
event <- anno_gen_1_2$event2
MN_info <- anno_gen_1_2$MN.Daughter
fam_vals <- unlist(rle(families)$values)
fam_reps <- as.numeric(rle(families)$lengths)*(length(unique(TPM_bychr$chr))-length(exclude_chrs))
fam_ids_vals <- unlist(rle(family_ids)$values)
fam_ids_reps <- as.numeric(rle(family_ids)$lengths)*(length(unique(TPM_bychr$chr))-length(exclude_chrs))

chr_num <- str_remove(unique(TPM_bychr$chr), pattern = "chr")
chr_num[chr_num == "X"] <- 23
chr_num[chr_num == "10a"] <- 10.1; chr_num[chr_num == "10b"] <- 10.2;
chr_num <- as.numeric(chr_num)
chr_num <- unique(TPM_bychr$chr)

#Prepare Visual Objects
visual.data <- data.table(TPM = as.numeric(unlist(flatten(TPM_bychr[,..names]))))
visual.data$VAR_A <- as.numeric(unlist(flatten(VAR_A[,..names])))
visual.data$VAR_B <- as.numeric(unlist(flatten(VAR_B[,..names])))
visual.data$chr <- rep(unique(chr_num), length = length(visual.data$TPM))
visual.data$cell <- rep(names, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$ids <- rep(ids, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$family <- rep(fam_vals, fam_reps)
visual.data$family_ids <- rep(fam_ids_vals, fam_ids_reps)
visual.data$relationship <- rep(relationship, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$generation <- rep(gen, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$rupt_time <- rep(rupt_time, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$event <- rep(event, each=length(unique(chr_num))-length(exclude_chrs))
visual.data$MN_info <- rep(MN_info, each=length(unique(chr_num))-length(exclude_chrs))

######################
### Visualize Data ###
######################

#add pvals to visual data object
visual.data$pval_loss <- as.numeric(unlist(flatten(pval_matrix_loss_bychr[,..names])))
visual.data$pval_control <- as.numeric(unlist(flatten(pval_matrix_control_bychr[,..names])))
visual.data$pval_gain <- as.numeric(unlist(flatten(pval_matrix_gain_bychr[,..names])))

#merge cell lineage data with possible events
setkey(visual.data, "chr")
possible_data2 <- visual.data[,c("chr", "relationship", "cell", "ids", "family", "family_ids", "event", "generation", "rupt_time", "MN_info", "TPM", "VAR_A", "VAR_B", "pval_loss", "pval_control", "pval_gain")]

#get state classifications
alpha <- .05 #normal alpha for trisomy and disomy as those distributions as samples size is 1
alpha_bonf <- .05/length(unique(visual.data$chr)) #Bonferroni corrected pval (due to 20 chrs)
loss_rows <- intersect( intersect( which(possible_data2$pval_loss >= alpha), which(possible_data2$pval_control < alpha_bonf)), which(possible_data2$pval_gain < alpha) )
loss_table <- cbind( rep("monosomy", length(loss_rows)), loss_rows)
control_rows <- intersect( intersect( which(possible_data2$pval_loss < alpha), which(possible_data2$pval_control >= alpha_bonf)), which(possible_data2$pval_gain < alpha) )
control_table <- cbind( rep("disomy", length(control_rows)), control_rows)
gain_rows <- intersect( intersect( which(possible_data2$pval_loss < alpha), which(possible_data2$pval_control < alpha_bonf)), which(possible_data2$pval_gain >= alpha) )
gain_table <- cbind( rep("trisomy", length(gain_rows)), gain_rows)
ploidy_table <- rbind(rbind(loss_table, control_table), gain_table)
ploidy_table <- data.table(ploidy_table[,c(2,1)]); colnames(ploidy_table) <- c("row", "state")
ploidy_table$row <- as.numeric(ploidy_table$row)
setkey(ploidy_table, "row")

#fill intermediate states
index <- data.table(row = (1:length(possible_data2$TPM)))
setkey(index, "row")
ploidy_table2 <- merge(index, ploidy_table, all.x = T) 
ploidy_table2$state[is.na(ploidy_table2$state)] <- "intermediate"

#combine state annotations
visual.data2 <- cbind(possible_data2, state = as.factor(ploidy_table2$state))
keycol <-c("generation", "family", "chr","relationship")
setorderv(visual.data2, keycol)

intermediate_32_rows <- intersect(which(visual.data2$state == "intermediate"), intersect(which(visual.data2$TPM > 1), which(visual.data2$TPM < 1.45) ) )
visual.data2[intermediate_32_rows,]$state <- "intermediate_32"
intermediate_12_rows <- intersect(which(visual.data2$state == "intermediate"), intersect(which(visual.data2$TPM < 1), which(visual.data2$TPM > .55) ) )
visual.data2[intermediate_12_rows,]$state <- "intermediate_12"
high_rows <- intersect(which(visual.data2$state == "intermediate"), which(visual.data2$TPM >= high_lim) )
visual.data2[high_rows,]$state <- "high"
low_rows <- intersect(which(visual.data2$state == "intermediate"), which(visual.data2$TPM <= low_lim) )
visual.data2[low_rows,]$state <- "low"

#Add a row for the MN chromosome
MN_fam_rle <- rle(visual.data2$family_ids)
visual.data2$MN_Cell <- as.character(rep(NA, nrow(visual.data2)))
visual.data2$MN_Sister <- as.character(rep(NA, nrow(visual.data2)))
for(i in 1:length(MN_fam_rle$values)) {
  cur_fam <- visual.data2[which(visual.data2$family_ids == MN_fam_rle$values[i]),]
  if(cur_fam$MN_info[1] == "A") {alt_info <- "B"}
  if(cur_fam$MN_info[1] == "B") {alt_info <- "A"}
  if(cur_fam$MN_info[1] == "AB") {alt_info <- "AB"}
  id <- unique(cur_fam$ids[which(cur_fam$relationship == cur_fam$MN_info)])
  id_alt <- unique(cur_fam$ids[which(cur_fam$relationship == alt_info)])
  if(length(id) != 0) {visual.data2$MN_Cell[which(visual.data2$family_ids == MN_fam_rle$values[i])] <- paste(unique(cur_fam$MN_info), ":", id)}
  if(length(id_alt) != 0) {visual.data2$MN_Sister[which(visual.data2$family_ids == MN_fam_rle$values[i])] <- paste(alt_info, ":", id_alt)}
}

#Save visual.data2 storing all classification data
saveRDS(visual.data2, file = sprintf("%s/aggregated_results/classification_data.rds", dirpath))


#cast the data into table
seg_table_TPM <- dcast(visual.data2, generation + family + family_ids + MN_Cell + MN_Sister + MN_info + event + rupt_time + chr ~ relationship, value.var = "TPM")[,c("A", "B", "c1", "c2")]
colnames(seg_table_TPM) <- c("A_TPM", "B_TPM", "c1_TPM", "c2_TPM") 
seg_table_VAR_A <- dcast(visual.data2, generation + family + family_ids + MN_Cell + MN_Sister + MN_info + event + rupt_time + chr ~ relationship, value.var = "VAR_A")[,c("A", "B", "c1", "c2")] #, id.vars = c("family", "chr"), measure.vars = c("state"))
colnames(seg_table_VAR_A) <- c("A_AlleleA", "B_AlleleA", "c1_AlleleA", "c2_AlleleA") 
seg_table_VAR_B <- dcast(visual.data2, generation + family + family_ids + MN_Cell + MN_Sister + MN_info + event + rupt_time + chr ~ relationship, value.var = "VAR_B")[,c("A", "B", "c1", "c2")] #, id.vars = c("family", "chr"), measure.vars = c("state"))
colnames(seg_table_VAR_B) <- c("A_AlleleB", "B_AlleleB", "c1_AlleleB", "c2_AlleleB") 
seg_table <- dcast(visual.data2, generation + family + family_ids + MN_Cell + MN_Sister + MN_info + event + rupt_time + chr ~ relationship, value.var = "state") #, id.vars = c("family", "chr"), measure.vars = c("state"))
seg_table <- cbind(seg_table[,c("generation","family","family_ids","MN_Cell","MN_Sister","MN_info","event","rupt_time","chr")], cbind(seg_table_TPM, cbind(seg_table_VAR_A, cbind(seg_table_VAR_B, seg_table[,c("A", "B", "c1", "c2")]) )))
write.csv(seg_table, file = sprintf("%s/aggregated_results/all_chr_profiles.csv", dirpath), row.names = F)

#############################################
### Get information for cells with events ###
#############################################

#remove cells that are completely diploid
cols <- (ncol(seg_table)-3):(ncol(seg_table))
rows_to_remove3 <- c()
for (i in 1:nrow(seg_table)) {if( length(setdiff(unique(seg_table[i,cols][!is.na(seg_table[i,cols])]), c("disomy"))) < 1 ){rows_to_remove3 <- append(rows_to_remove3, i)}}
seg_table_allevents <- seg_table[-rows_to_remove3,]
seg_table_diploid_events <- seg_table[rows_to_remove3,]
write.csv(seg_table_allevents, file = sprintf("%s/aggregated_results/all_events.csv", dirpath), row.names = F)
write.csv(seg_table_diploid_events, file = sprintf("%s/aggregated_results/diploid_events.csv", dirpath), row.names = F)

#Add back single row for families that are completely diploid
rows_to_keep <- c()
for (fam in unique(seg_table$family)) {
  cur_fam <- seg_table[which(seg_table$family == fam),]
  cntr <- 0
  for (i in 1:nrow(cur_fam)) {
    if ( length(setdiff(cur_fam[i,cols][!is.na(cur_fam[i,cols])], c("disomy"))) > 0 ) {cntr <- cntr + 1}
  }
  if (cntr == 0) {rows_to_keep <- append(rows_to_keep, which(seg_table$family == fam))}
}
seg_table_diploid_families <- seg_table[rows_to_keep,]
seg_table_diploid_families_red <- seg_table_diploid_families[which(seg_table_diploid_families$chr == "chr1"),]
seg_table_diploid_families_red$chr <- "all chrs"
write.csv(seg_table_diploid_families_red, file = sprintf("%s/aggregated_results/diploid_families.csv", dirpath), row.names = F)

#add all events plus place holder row for completely diploid families
all_family_events_2 <- data.table(rbind(seg_table_allevents, seg_table_diploid_families_red))

#reorder columns for easier analysis
all_family_events_2$chr[all_family_events_2$chr == "all chrs"] <- NA; all_family_events_2$chr[all_family_events_2$chr == "chrX"] <- "chr23";
all_family_events_2$chr[all_family_events_2$chr == "chr10a"] <- "chr10.1"; all_family_events_2$chr[all_family_events_2$chr == "chr10b"] <- "chr10.2";
all_family_events_2$chr <- as.numeric(str_replace(string = all_family_events_2$chr, pattern = "chr", replacement = ""))
setkey(all_family_events_2, "generation", "family", "chr")

########################################
### Multiply TPM ratios columns by 2 ###
########################################

#multiply these 4 rows by 2
changeCols <- c("A_TPM", "B_TPM", "c1_TPM", "c2_TPM") 
all_family_events_2[,(changeCols):= lapply(.SD, function(x){x*2}), .SDcols = changeCols] 

#############################
### Annotate MN haplotype ###
#############################

#if else representation of MN hap reasoning
hap_column <- c()
for(i in 1:nrow(all_family_events_2)) {
  tmp_hap <- c()
  tmp_event <- all_family_events_2[i,]
  if(tmp_event$chr %in% c(10.2, 23)){
    tmp_event_alt <- all_family_events_2[which(all_family_events_2$family == tmp_event$family),]
    tmp_event_alt_2 <- tmp_event_alt[which(tmp_event_alt$chr == 10.2),]
    if(nrow(tmp_event_alt_2) > 0){
      if(tmp_event_alt_2$A %in% c("monosomy", "intermediate_12", "low") | tmp_event$B %in% c("monosomy", "intermediate_12", "low")) {
        if(tmp_event_alt_2$A %in% c("monosomy", "intermediate_12", "low")){
          if(tmp_event_alt_2$A_AlleleA < (tmp_event_alt_2$A_AlleleB - .33)){tmp_hap <- append(tmp_hap, "A")}
          if(tmp_event_alt_2$A_AlleleA > (tmp_event_alt_2$A_AlleleB - .33)){tmp_hap <- append(tmp_hap, "B")}
        }
        if(tmp_event_alt_2$B %in% c("monosomy", "intermediate_12", "low")){
          if(tmp_event_alt_2$B_AlleleA < (tmp_event_alt_2$B_AlleleB - .33)){tmp_hap <- append(tmp_hap, "A")}
          if(tmp_event_alt_2$B_AlleleA > (tmp_event_alt_2$B_AlleleB - .33)){tmp_hap <- append(tmp_hap, "B")}
        }
      }
      if(tmp_event_alt_2$c1 %in% c("monosomy", "intermediate_12", "low") & tmp_event_alt_2$c2 %in% c("monosomy", "intermediate_12", "low")){
        if(tmp_event_alt_2$c1_AlleleA < (tmp_event_alt_2$c1_AlleleB - .33)){tmp_hap <- append(tmp_hap, "A")}
        if(tmp_event_alt_2$c1_AlleleA > (tmp_event_alt_2$c1_AlleleB - .33)){tmp_hap <- append(tmp_hap, "B")}
      }
      #I did not account for 3:2:3:3 cases in X or 10b (doesn't occurr and complicated implementation)
    }
  } else {
    if(tmp_event$A %in% c("monosomy", "intermediate_12", "low") | tmp_event$B %in% c("monosomy", "intermediate_12", "low")) {
      if(tmp_event$A %in% c("monosomy", "intermediate_12", "low")){
        if(tmp_event$A_AlleleA < tmp_event$A_AlleleB){tmp_hap <- append(tmp_hap, "A")}
        if(tmp_event$A_AlleleA > tmp_event$A_AlleleB){tmp_hap <- append(tmp_hap, "B")}
      }
      if(tmp_event$B %in% c("monosomy", "intermediate_12", "low")){
        if(tmp_event$B_AlleleA < tmp_event$B_AlleleB){tmp_hap <- append(tmp_hap, "A")}
        if(tmp_event$B_AlleleA > tmp_event$B_AlleleB){tmp_hap <- append(tmp_hap, "B")}
      }
    }
    if(tmp_event$c1 %in% c("monosomy", "intermediate_12", "low") & tmp_event$c2 %in% c("monosomy", "intermediate_12", "low")){
      if(tmp_event$c1_AlleleA < tmp_event$c1_AlleleB){tmp_hap <- append(tmp_hap, "A")}
      if(tmp_event$c1_AlleleA > tmp_event$c1_AlleleB){tmp_hap <- append(tmp_hap, "B")}
    }
    if(tmp_event$c1 %in% c("trisomy", "intermediate_32", "high") & tmp_event$c2 %in% c("trisomy", "intermediate_32", "high")){
      if(tmp_event$c1_AlleleA > tmp_event$c1_AlleleB){tmp_hap <- append(tmp_hap, "A")}
      if(tmp_event$c1_AlleleA < tmp_event$c1_AlleleB){tmp_hap <- append(tmp_hap, "B")}
    }
    if(tmp_event$generation == "Gen1" & tmp_event$A %in% c("disomy") & tmp_event$B %in% c("trisomy", "intermediate_32", "high")){
      if(tmp_event$B_AlleleA > tmp_event$B_AlleleB){tmp_hap <- append(tmp_hap, "A")}
      if(tmp_event$B_AlleleA < tmp_event$B_AlleleB){tmp_hap <- append(tmp_hap, "B")}
    }
    if(tmp_event$generation == "Gen2" & tmp_event$A %in% c("trisomy", "intermediate_32", "high") & tmp_event$B %in% c("disomy") & is.na(tmp_event$c1) & is.na(tmp_event$c2)){
      if(tmp_event$A_AlleleA > tmp_event$A_AlleleB){tmp_hap <- append(tmp_hap, "A")}
      if(tmp_event$A_AlleleA < tmp_event$A_AlleleB){tmp_hap <- append(tmp_hap, "B")}
    }
  }
  if(length(tmp_hap) > 1) {
    if(tmp_hap[1] == tmp_hap[2]) {tmp_hap <- tmp_hap[1]} else {print("Could not determine MN haplotpye"); tmp_hap <- NA;}
  }
  if(is.null(tmp_hap)) {hap_column <- append(hap_column, NA)} else {hap_column <- append(hap_column, tmp_hap)}
}

#add hap annotations to df
all_family_events_2$MN_hap <- hap_column

#############################################################
### add annotatation for which chrs are used as reference ###
#############################################################

fam_column <- c()
for(i in 1:nrow(hand_samples)){fam_column <- append(fam_column, anno$Pairs[which(anno$WTA.plate == hand_samples$ID[i])])}
hand_samples$fam <- fam_column

ref_column <- c()
for(j in 1:nrow(all_family_events_2)){
  tmp <- all_family_events_2[j,]
  tmp_mat <- hand_samples[which(hand_samples$fam == tmp$family),]
  tmp_cntr <- 0
  if(nrow(tmp_mat) > 0) {
    for(i in 1:nrow(tmp_mat)){
      if(tmp_mat$chr[i] == tmp$chr){tmp_cntr <- tmp_cntr + 1}
    }
  }
  if(tmp_cntr > 0){ref_column <- append(ref_column, TRUE)} else {ref_column <- append(ref_column, FALSE)}
}

all_family_events_2$reference <- ref_column

######################################
### Estimate CN pattern annotation ###
######################################

CN_values <- round(all_family_events_2[,c("A_TPM", "B_TPM", "c1_TPM", "c2_TPM")])
CN_pattern_column <- paste(CN_values$A_TPM, CN_values$B_TPM, CN_values$c1_TPM, CN_values$c2_TPM)
CN_pattern_column_2 <- sub(sub(sub(x = CN_pattern_column, pattern = "NA", replacement = ""), pattern = "NA", replacement = ""), pattern = "NA", replacement = "")

all_family_events_2$CN_pattern <- CN_pattern_column_2

###############################
### Hap specific CN pattern ###
###############################

hap_spec_CN_pattern_column <- data.table()
for(i in 1:nrow(all_family_events_2)){
  if(!is.na(all_family_events_2$MN_hap[i])){
    haps <- c(sprintf("A_Allele%s", all_family_events_2$MN_hap[i]), sprintf("B_Allele%s", all_family_events_2$MN_hap[i]), sprintf("c1_Allele%s", all_family_events_2$MN_hap[i]), sprintf("c2_Allele%s", all_family_events_2$MN_hap[i]))
    hap_spec_CN <- round(all_family_events_2[i,..haps])
    colnames(hap_spec_CN) <- c("A", "B", "c1", "c2")
    hap_spec_CN_pattern_column <- rbind(hap_spec_CN_pattern_column, hap_spec_CN)
  } else {hap_spec_CN_pattern_column <-  rbind(hap_spec_CN_pattern_column, data.table(NA, NA, NA, NA), use.names=FALSE)}
}
hap_spec_CN_pattern <- paste(hap_spec_CN_pattern_column$A, hap_spec_CN_pattern_column$B, hap_spec_CN_pattern_column$c1, hap_spec_CN_pattern_column$c2)
hap_spec_CN_pattern_2 <- sub(sub(sub(sub(x = hap_spec_CN_pattern, pattern = "NA", replacement = ""), pattern = "NA", replacement = ""), pattern = "NA", replacement = ""), pattern = "NA", replacement = "")

all_family_events_2$MNhap_CN_pattern <- hap_spec_CN_pattern_2

##########################
### AS_TPM for MN cell ###
##########################
hap_spec_CN_column <- c()
for(i in 1:nrow(all_family_events_2)){
  if(!is.na(all_family_events_2$MN_hap[i])){
    hap <- c(sprintf("A_Allele%s", all_family_events_2$MN_hap[i]))
    hap_spec_CN_column <- append(hap_spec_CN_column, all_family_events_2[i,..hap])
  } else {hap_spec_CN_column <- append(hap_spec_CN_column, NA)}
}

all_family_events_2$MNhap_TPM <- as.character(hap_spec_CN_column)

###############################
### Calculate hap CN state  ###
###############################

CN_state_column <- c()
for (i in 1:nrow(all_family_events_2)) {
  CN_state <- NA
  hap <- all_family_events_2$MN_hap[i]
  cell_id <- strsplit(all_family_events_2$MN_Cell[i], split= " : ")[[1]][2]
  cell <- anno$WTA.plate[which(anno$Cell_IDs == cell_id)]
  chr <- paste0("chr", all_family_events_2$chr[i])
  if(chr == "chr23"){chr <- "chrX"}; if(chr == "chr10.1"){chr <- "chr10a"}; if(chr == "chr10.2"){chr <- "chr10b"};
  MNhap_TPM <- as.numeric(all_family_events_2$MNhap_TPM[i])
  if(!chr %in% c("chrX", "chr10b")){
    if(!is.na(hap)) {
      if(hap == "A") {
        pvals_list <- c(AS_TPM_bychr_mono_pvals_chr_spec$A[,..cell][which(AS_TPM_bychr_mono_pvals_chr_spec$A$chr == chr)], AS_TPM_bychr_ctrl_pvals_chr_spec$A[,..cell][which(AS_TPM_bychr_ctrl_pvals_chr_spec$A$chr == chr)])
        CN_state <- c("normal", "gain")[which(pvals_list > alpha)]
        if(!is.na(MNhap_TPM)){
          if(MNhap_TPM > high_lim_AS_TPM){CN_state <- "high"}
          if(MNhap_TPM < low_lim_AS_TPM){CN_state <- "low"}
        }
        if(is_empty(CN_state)){if(MNhap_TPM<1){CN_state <- "loss"}else{CN_state <- "intermediate"}}
      } else {
        pvals_list <- c(AS_TPM_bychr_mono_pvals_chr_spec$B[,..cell][which(AS_TPM_bychr_mono_pvals_chr_spec$B$chr == chr)], AS_TPM_bychr_ctrl_pvals_chr_spec$B[,..cell][which(AS_TPM_bychr_ctrl_pvals_chr_spec$B$chr == chr)])
        CN_state <- c("normal", "gain")[which(pvals_list > alpha)]
        if(!is.na(MNhap_TPM)){
          if(MNhap_TPM > high_lim_AS_TPM){CN_state <- "high"}
          if(MNhap_TPM < low_lim_AS_TPM){CN_state <- "low"}
        }
        if(is_empty(CN_state)){if(MNhap_TPM<1){CN_state <- "loss"}else{CN_state <- "intermediate"}}
      }
    }
  }
  if(length(CN_state) > 1){print(CN_state)}
  if(length(CN_state) == 1){CN_state_column <- append(CN_state_column, CN_state)} else {CN_state_column <- append(CN_state_column, NA); break}
}

all_family_events_2$CN_state <- CN_state_column

##########################################
### Pull monosomy pvals for MNhap_TPMs ###
##########################################

chr_spec_mono_pval_column <- c()
for (i in 1:nrow(all_family_events_2)) {
  pvals_list <- c()
  hap <- all_family_events_2$MN_hap[i]
  cell_id <- strsplit(all_family_events_2$MN_Cell[i], split= " : ")[[1]][2]
  cell <- anno$WTA.plate[which(anno$Cell_IDs == cell_id)]
  chr <- paste0("chr", all_family_events_2$chr[i])
  if(chr == "chr23"){chr <- "chrX"}; if(chr == "chr10.1"){chr <- "chr10a"}; if(chr == "chr10.2"){chr <- "chr10b"};
  MNhap_TPM <- all_family_events_2$MNhap_TPM[i]
  if(!is.na(hap)) {
    if(hap == "A") {
      pvals_list <- c(AS_TPM_bychr_mono_pvals_chr_spec$A[,..cell][which(AS_TPM_bychr_mono_pvals_chr_spec$A$chr == chr)])
    } else {
      pvals_list <- c(AS_TPM_bychr_mono_pvals_chr_spec$B[,..cell][which(AS_TPM_bychr_mono_pvals_chr_spec$B$chr == chr)])
    }
  }
  pval <- unname(unlist(pvals_list))
  if(length(pval) > 1){print(i)}
  if(length(pval) == 1){chr_spec_mono_pval_column <- append(chr_spec_mono_pval_column, pval)} else {chr_spec_mono_pval_column <- append(chr_spec_mono_pval_column, NA)}
}

all_family_events_2$chr_spec_mono_pval <- chr_spec_mono_pval_column

########################################
### Pull disomy pvals for MNhap_TPMs ###
########################################

chr_spec_ctrl_pval_column <- c()
for (i in 1:nrow(all_family_events_2)) {
  pvals_list <- c()
  hap <- all_family_events_2$MN_hap[i]
  cell_id <- strsplit(all_family_events_2$MN_Cell[i], split= " : ")[[1]][2]
  cell <- anno$WTA.plate[which(anno$Cell_IDs == cell_id)]
  chr <- paste0("chr", all_family_events_2$chr[i])
  if(chr == "chr23"){chr <- "chrX"}; if(chr == "chr10.1"){chr <- "chr10a"}; if(chr == "chr10.2"){chr <- "chr10b"};
  MNhap_TPM <- all_family_events_2$MNhap_TPM[i]
  if(!is.na(hap)) {
    if(hap == "A") {
      pvals_list <- c(AS_TPM_bychr_ctrl_pvals_chr_spec$A[,..cell][which(AS_TPM_bychr_ctrl_pvals_chr_spec$A$chr == chr)])
    } else {
      pvals_list <- c(AS_TPM_bychr_ctrl_pvals_chr_spec$B[,..cell][which(AS_TPM_bychr_ctrl_pvals_chr_spec$B$chr == chr)])
    }
  }
  pval <- unname(unlist(pvals_list))
  if(length(pval) > 1){print(i)}
  if(length(pval) == 1){chr_spec_ctrl_pval_column <- append(chr_spec_ctrl_pval_column, pval)} else {chr_spec_ctrl_pval_column <- append(chr_spec_ctrl_pval_column, NA)}
}

all_family_events_2$chr_spec_ctrl_pval <- chr_spec_ctrl_pval_column

##########################################
### Reorder, round and save all events ###
##########################################

#Reorder cols
ordered_cols <- c("generation","family","family_ids","MN_Cell","MN_Sister","MN_info","event","rupt_time","chr","A_TPM","A","A_AlleleA","A_AlleleB","B_TPM","B","B_AlleleA","B_AlleleB","c1_TPM","c1","c1_AlleleA","c1_AlleleB","c2_TPM","c2","c2_AlleleA","c2_AlleleB","CN_pattern", "MN_hap", "MNhap_CN_pattern", "MNhap_TPM", "chr_spec_mono_pval", "chr_spec_ctrl_pval", "CN_state", "reference")
all_family_events_2 <- all_family_events_2[,..ordered_cols]
setkey(all_family_events_2, "generation", "family", "chr")

#save event table
write.csv(all_family_events_2, file = sprintf("%s/aggregated_results/all_family_events.csv", dirpath), row.names = F)
all_family_events_2 <- read.csv(file = sprintf("%s/aggregated_results/all_family_events.csv", dirpath))
