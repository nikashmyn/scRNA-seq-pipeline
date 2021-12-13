#automated analysis script

########################
### Passed Arguments ###
########################

#TODO: check if there are any completely diploid families in the analysis and then
#          add code to remove the "rows_to_keep" from the "rows_to_remove"
#some intermediate values still in visuals results
#double intermediate_32 or intermediate_21 in gen1 should be removed. 

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
print(args)
scriptsdir <- args[1]
wkpath <- dirpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", wkpath)
system(mk_agg_dir)

####################
### Read-in data ###
####################

library(reshape)

#Read in data
golden_samples <- readRDS(file = sprintf("%s/aggregated_results/golden_set_tpms_bychr.rds", dirpath))

#Read in grouped sample information
hand_samples <- readRDS(sprintf("%s/aggregated_results/grouped_control_aneuploidies.rds", dirpath))

#read in normed TPM by chr
adt.bychr <- readRDS(file=sprintf("%s/aggregated_results/adt.bychr.rds", dirpath))

#Allele fraction by chr
allele_frac.bychr <- readRDS(file=sprintf("%s/aggregated_results/allele_frac.bychr.rds", dirpath))
VAR_A <- copy(allele_frac.bychr$A)
VAR_B <- copy(allele_frac.bychr$B)

#read in pvalues from control dist bychr
pval_matrix_control_bychr <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_control_bychr.rds", dirpath))
pval_matrix_gain_bychr <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_gain_bychr.rds", dirpath))
pval_matrix_loss_bychr <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_loss_bychr.rds", dirpath))

#Read in annotation list with sister and cousin information
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))

####################
### Boxplot Data ###
####################

#main data to be plotted byarm
boxplots.vals <- cbind(golden_samples$loss_bychr$tpm, rep("loss", length(golden_samples$loss_bychr$tpm))); 
boxplots.vals <- rbind(boxplots.vals, cbind(golden_samples$control_bychr$tpm, rep("control", length(golden_samples$control_bychr$tpm))))
boxplots.vals <- data.table(rbind(boxplots.vals, cbind(golden_samples$gain_bychr$tpm, rep("gain", length(golden_samples$gain_bychr$tpm)))))
colnames(boxplots.vals) <- c("TPM", "Group")
boxplots.vals$TPM <- as.numeric(boxplots.vals$TPM)
boxplots.vals$Group <- as.factor(boxplots.vals$Group)
dist_colors <- c(rep(1, length(golden_samples$loss_bychr$tpm)), rep(10, length(golden_samples$control_bychr$tpm)), rep(20, length(golden_samples$gain_bychr$tpm)))
boxplots.vals$color <- c(dist_colors)

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
fam_reps <- as.numeric(rle(families)$lengths)*(length(unique(adt.bychr$seqnames))-length(exclude_chrs))
fam_ids_vals <- unlist(rle(family_ids)$values)
fam_ids_reps <- as.numeric(rle(family_ids)$lengths)*(length(unique(adt.bychr$seqnames))-length(exclude_chrs))

#Prepare Visual Objects
visual.data <- data.table(TPM = as.numeric(flatten(adt.bychr[,..names])))
visual.data$VAR_A <- as.numeric(flatten(VAR_A[,..names]))
visual.data$VAR_B <- as.numeric(flatten(VAR_B[,..names]))
visual.data$chr <- as.numeric(rep(unique(adt.bychr$seqnames), length = length(visual.data$TPM)))
visual.data$cell <- rep(names, each=length(unique(adt.bychr$seqnames))-length(exclude_chrs))
visual.data$ids <- rep(ids, each=length(unique(adt.bychr$seqnames))-length(exclude_chrs))
visual.data$family <- rep(fam_vals, fam_reps)
visual.data$family_ids <- rep(fam_ids_vals, fam_ids_reps)
visual.data$relationship <- rep(relationship, each=length(unique(adt.bychr$seqnames))-length(exclude_chrs))
visual.data$generation <- rep(gen, each=length(unique(adt.bychr$seqnames))-length(exclude_chrs))
visual.data$rupt_time <- rep(rupt_time, each=length(unique(adt.bychr$seqnames))-length(exclude_chrs))
visual.data$event <- rep(event, each=length(unique(adt.bychr$seqnames))-length(exclude_chrs))
visual.data$MN_info <- rep(MN_info, each=length(unique(adt.bychr$seqnames))-length(exclude_chrs))

######################
### Visualize Data ###
######################

#Plot the boxplots and strip plots
boxplot.visual <- ggplot() + 
  geom_jitter(data = visual.data, mapping = aes(x = 4, y = TPM, color = chr), width = 0.3) +
  geom_boxplot(data = boxplots.vals, mapping = aes(x = Group, y = TPM, color = color), outlier.alpha = 0)  +
  ylim(c(.25,1.75)) + scale_x_discrete(limits = c("1", "2", "3", "4"), labels=c("loss", "control", "gain", "family")) +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Classes | %s ", myfamily))

#add pvals to visual data object
visual.data$pval_loss <- as.numeric(flatten(pval_matrix_loss_bychr[,..names]))
visual.data$pval_control <- as.numeric(flatten(pval_matrix_control_bychr[,..names]))
visual.data$pval_gain <- as.numeric(flatten(pval_matrix_gain_bychr[,..names]))

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
high_rows <- intersect(which(visual.data2$state == "intermediate"), which(visual.data2$TPM >= 1.45) )
visual.data2[high_rows,]$state <- "high"
low_rows <- intersect(which(visual.data2$state == "intermediate"), which(visual.data2$TPM <= .55) )
visual.data2[low_rows,]$state <- "low"

#Add a row for the MN chromosome
MN_fam_rle <- rle(visual.data2$family_ids)
visual.data2$MN_Cell <- as.character(rep(NA, nrow(visual.data2)))
visual.data2$MN_Sister <- as.character(rep(NA, nrow(visual.data2)))
for(i in 1:length(MN_fam_rle$values)) {
  cur_fam <- visual.data2[which(visual.data2$family_ids == MN_fam_rle$values[i]),]
  if(cur_fam$MN_info == "A") {alt_info <- "B"}
  if(cur_fam$MN_info == "B") {alt_info <- "A"}
  if(cur_fam$MN_info == "AB") {alt_info <- "AB"}
  id <- unique(cur_fam$ids[which(cur_fam$relationship == cur_fam$MN_info)])
  id_alt <- unique(cur_fam$ids[which(cur_fam$relationship == alt_info)])
  if(length(id) != 0) {visual.data2$MN_Cell[which(visual.data2$family_ids == MN_fam_rle$values[i])] <- paste(unique(cur_fam$MN_info), ":", id)}
  if(length(id_alt) != 0) {visual.data2$MN_Sister[which(visual.data2$family_ids == MN_fam_rle$values[i])] <- paste(alt_info, ":", id_alt)}
}

#cast the data into table
seg_table_TPM <- dcast(visual.data2, generation + family + family_ids + MN_Cell + MN_Sister + MN_info + event + rupt_time + chr ~ relationship, value.var = "TPM")[,c("A", "B", "c1", "c2")]
colnames(seg_table_TPM) <- c("A_TPM", "B_TPM", "c1_TPM", "c2_TPM") 
seg_table_VAR_A <- dcast(visual.data2, generation + family + family_ids + MN_Cell + MN_Sister + MN_info + event + rupt_time + chr ~ relationship, value.var = "VAR_A")[,c("A", "B", "c1", "c2")] #, id.vars = c("family", "chr"), measure.vars = c("state"))
colnames(seg_table_VAR_A) <- c("A_AlleleA", "B_AlleleA", "c1_AlleleA", "c2_AlleleA") 
seg_table_VAR_B <- dcast(visual.data2, generation + family + family_ids + MN_Cell + MN_Sister + MN_info + event + rupt_time + chr ~ relationship, value.var = "VAR_B")[,c("A", "B", "c1", "c2")] #, id.vars = c("family", "chr"), measure.vars = c("state"))
colnames(seg_table_VAR_B) <- c("A_AlleleB", "B_AlleleB", "c1_AlleleB", "c2_AlleleB") 
seg_table <- dcast(visual.data2, generation + family + family_ids + MN_Cell + MN_Sister + MN_info + event + rupt_time + chr ~ relationship, value.var = "state") #, id.vars = c("family", "chr"), measure.vars = c("state"))
seg_table <- cbind(seg_table[,c("generation","family","family_ids","MN_Cell","MN_Sister","MN_info","event","rupt_time","chr")], cbind(seg_table_TPM, cbind(seg_table_VAR_A, cbind(seg_table_VAR_B, seg_table[,c("A", "B", "c1", "c2")]) )))
write_csv(seg_table, file = sprintf("%s/aggregated_results/all_chr_profiles.csv", dirpath))

###################################
### semi-automate gen1 analysis ###
###################################

#let just look at gen1
seg_table_gen1 <- seg_table[which(seg_table$generation == "Gen1"),]

#get cells that are completely diploid
cols <- (ncol(seg_table_gen1)-3):(ncol(seg_table_gen1))
rows_to_keep <- c()
for (fam in unique(seg_table_gen1$family)) {
  cur_fam <- seg_table_gen1[which(seg_table_gen1$family == fam),]
  cntr <- 0
  for (i in 1:nrow(cur_fam)) {
    if ( length(setdiff(cur_fam[i,cols][!is.na(cur_fam[i,cols])], c("disomy"))) > 0 ) {cntr <- cntr + 1}
  }
  if (cntr == 0) {rows_to_keep <- append(rows_to_keep, which(seg_table_gen2$family == fam))}
}

#get rid of completely diploid chrs
scenarios <- seg_table_gen1[,cols]; rows_to_remove1 <- c()
for (i in 1:nrow(seg_table_gen1)) {if( length(setdiff(unique(scenarios[i,][!is.na(scenarios[i,])]), c("disomy"))) < 1 ){rows_to_remove1 <- append(rows_to_remove1, i)}}

#remove events where all cells are mono. In the case of 2 cell fams you could only get a double mono with a persistant MN which we do not have in the two cell families. 
#for (i in 1:nrow(seg_table_gen2)) {if( length(setdiff(unique(scenarios[i,][!is.na(scenarios[i,])]), c("monosomy"))) < 1 ){rows_to_remove1 <- append(rows_to_remove1, i)}}
#for (i in 1:nrow(seg_table_gen2)) {if( length(setdiff(unique(scenarios[i,][!is.na(scenarios[i,])]), c("monosomy"))) < 1 ){print(scenarios[i,])}}

#remove events where there are >= 2 trisomies in a family with less than 5 cells as these cannot be interpretable cases
for (i in 1:nrow(seg_table_gen1)) { 
  if(length(scenarios[i,][!is.na(scenarios[i,])]) < 5) {
    if( length(which(scenarios[i,] == "trisomy")) >= 2 )  {
      rows_to_remove1 <- append(rows_to_remove1, i)}}}

#remove events where there are only disomies and trisomies
for (i in 1:nrow(seg_table_gen1)) { 
  if( sum(unique(scenarios[i,][!is.na(scenarios[i,])]) %in% c("disomy", "intermediate_32", "trisomy")) >= 2 )  {
    rows_to_remove1 <- append(rows_to_remove1, i)}}

#remove events where 2:32, 32:32, 21:21
for (i in 1:nrow(seg_table_gen1)) { 
  if( length(setdiff(scenarios[i,][!is.na(scenarios[i,])], c("disomy", "intermediate_32"))) < 1 ) {
    rows_to_remove1 <- append(rows_to_remove1, i)}
  if( length(setdiff(scenarios[i,][!is.na(scenarios[i,])], c("intermediate_12"))) < 1 ) {
    rows_to_remove1 <- append(rows_to_remove1, i)}
}

#get remaining possible events after event criteria
rows_to_remove1 <- unique(rows_to_remove1)
seg_table_gen1_reduced <- seg_table_gen1[-rows_to_remove1,]

#Finalize samples with only one event
one_event_families <- rle(seg_table_gen1_reduced$family)$values[which(rle(seg_table_gen1_reduced$family)$lengths <= 1)]
seg_table_gen1_final <- seg_table_gen1_reduced[which(seg_table_gen1_reduced$family %in% one_event_families),]

#get the families with more than one possible event
seg_table_gen1_mult_event <- seg_table_gen1_reduced[which(!seg_table_gen1_reduced$family %in% one_event_families),]

#save results as csv
write_csv(seg_table_gen1_reduced, file = sprintf("%s/aggregated_results/possible_MN_events_gen1.csv", dirpath))



###################################
### semi-automate gen2 analysis ###
###################################

#let just look at gen2
seg_table_gen2 <- seg_table[which(seg_table$generation == "Gen2"),]

#get cells that are completely diploid
rows_to_keep <- c()
for (fam in unique(seg_table_gen2$family)) {
  cur_fam <- seg_table_gen2[which(seg_table_gen2$family == fam),]
  cntr <- 0
  for (i in 1:nrow(cur_fam)) {
    if ( length(setdiff(cur_fam[i,cols][!is.na(cur_fam[i,cols])], c("disomy"))) > 0 ) {cntr <- cntr + 1}
  }
  if (cntr == 0) {rows_to_keep <- append(rows_to_keep, which(seg_table_gen2$family == fam))}
}

#get rid of completely diploid chrs
scenarios <- seg_table_gen2[,cols]; rows_to_remove2 <- c()
for (i in 1:nrow(seg_table_gen2)) {if( length(setdiff(unique(scenarios[i,][!is.na(scenarios[i,])]), c("disomy"))) < 1 ){rows_to_remove2 <- append(rows_to_remove2, i)}}

#remove events where all cells are mono. In the case of 2 cell fams you could only get a double mono with a persistant MN which we do not have in the two cell families. 
#for (i in 1:nrow(seg_table_gen2)) {if( length(setdiff(unique(scenarios[i,][!is.na(scenarios[i,])]), c("monosomy"))) < 1 ){rows_to_remove2 <- append(rows_to_remove2, i)}}
#for (i in 1:nrow(seg_table_gen2)) {if( length(setdiff(unique(scenarios[i,][!is.na(scenarios[i,])]), c("monosomy"))) < 1 ){print(scenarios[i,])}}

#remove events where cousins/internal controls are not the same ploidy. These are most likely noise events/uninterpretable.
for (i in 1:nrow(seg_table_gen2)) {
  if( length(unique(scenarios[i,c(3:4)][!is.na(scenarios[i,c(3:4)])])) > 1 ){
    if("disomy" %in% scenarios[i,c(3:4)][!is.na(scenarios[i,c(3:4)])]) {
      rows_to_remove2 <- append(rows_to_remove2, i)}}}

#remove events where there are >= 2 trisomies in a family with less than 5 cells as these cannot be interpretable cases
for (i in 1:nrow(seg_table_gen2)) { 
  if(length(scenarios[i,][!is.na(scenarios[i,])]) < 5) {
    if( length(which(scenarios[i,] == "trisomy")) >= 2 )  {
      rows_to_remove2 <- append(rows_to_remove2, i)}}}

#remove events where there are only disomies, trisomies, and intermediate_32s in a 3 cell fam
for (i in 1:nrow(seg_table_gen2)) { 
  if(length(scenarios[i,][!is.na(scenarios[i,])]) == 3) {
    if(length( setdiff(unique(scenarios[i,][!is.na(scenarios[i,])]), c("disomy", "intermediate_32", "trisomy"))) < 1) { 
      rows_to_remove2 <- append(rows_to_remove2, i)}}}

#remove events where cousins are disomic and sisters intermediate_32 and disomic
for (i in 1:nrow(seg_table_gen2)) { 
  if(length(scenarios[i,][!is.na(scenarios[i,])]) > 2) {
    if( length(setdiff(c("disomy"), scenarios[i,][!is.na(scenarios[i,])][-c(1:2)])) < 1 ) {
      if( length(setdiff(scenarios[i,][!is.na(scenarios[i,])][c(1:2)], c("disomy", "intermediate_32"))) < 1 ) {
        rows_to_remove2 <- append(rows_to_remove2, i)}}}}

#get remaining possible events after event criteria
rows_to_remove2 <- unique(rows_to_remove2)
seg_table_gen2_reduced <- seg_table_gen2[-rows_to_remove2,]

#Finalize samples with only one event
one_event_families <- rle(seg_table_gen2_reduced$family)$values[which(rle(seg_table_gen2_reduced$family)$lengths <= 1)]
seg_table_gen2_final <- seg_table_gen2_reduced[which(seg_table_gen2_reduced$family %in% one_event_families),]

#get the families with more than one possible event
seg_table_gen2_mult_event <- seg_table_gen2_reduced[which(!seg_table_gen2_reduced$family %in% one_event_families),]

#save results as csv
write_csv(seg_table_gen2_reduced, file = sprintf("%s/aggregated_results/possible_MN_events_gen2.csv", dirpath))

###########################################
### bind gen1 and gen2 results together ###
###########################################

seg_table_reduced <- rbind(seg_table_gen1_reduced, seg_table_gen2_reduced)
write_csv(seg_table_reduced, file = sprintf("%s/aggregated_results/possible_MN_events.csv", dirpath))

############################################
### segment chrs into MN and non MN chrs ###
############################################

#transform indices so they are correct to the both gen table
rows_to_remove2 <- rows_to_remove2 + nrow(seg_table[which(seg_table$generation == "Gen1"),])
rows_to_remove <- append(rows_to_remove1, rows_to_remove2)

#Get the compliment of the possible MN chrs
seg_table_noMN <- seg_table[rows_to_remove,]

################################################
### Get information for cells with no events ###
################################################

seg_table_noevent <- seg_table[which(!seg_table$family %in% unique(seg_table_reduced$family)),]
write_csv(seg_table_noevent, file = sprintf("%s/aggregated_results/cells_noMN_events.csv", dirpath))

################################################
### Get information for cells with no events ###
################################################

#remove cells that are completely diploid
rows_to_remove3 <- c()
for (i in 1:nrow(seg_table)) {if( length(setdiff(unique(seg_table[i,cols][!is.na(seg_table[i,cols])]), c("disomy"))) < 1 ){rows_to_remove3 <- append(rows_to_remove3, i)}}
seg_table_allevents <- seg_table[-rows_to_remove3,]
seg_table_diploid_events <- seg_table[rows_to_remove3,]
write_csv(seg_table_allevents, file = sprintf("%s/aggregated_results/all_events.csv", dirpath))
write_csv(seg_table_diploid_events, file = sprintf("%s/aggregated_results/diploid_events.csv", dirpath))

#get cells that are completely diploid
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
seg_table_diploid_families_red <- seg_table_diploid_families[which(seg_table_diploid_families$chr == 1),]
seg_table_diploid_families_red$chr <- "all chrs"
write_csv(seg_table_diploid_families_red, file = sprintf("%s/aggregated_results/diploid_families.csv", dirpath))

all_family_events <- rbind(seg_table_allevents, seg_table_diploid_families_red)
write_csv(all_family_events, file = sprintf("%s/aggregated_results/all_family_events.csv", dirpath))

#################################################
### get all data for sisters in results table ###
#################################################

rows_to_retrieve_compliment <- c()
for(i in 1:nrow(seg_table_noMN)) {
  chr <- seg_table_noMN[i,]$chr
  family <- seg_table_noMN[i,]$family
  index <- intersect(which(visual.data2$chr == chr), which(visual.data2$family == family))
  rows_to_retrieve_compliment <- append(rows_to_retrieve_compliment, index)
}

rows_to_retrieve <- c()
for(i in 1:nrow(seg_table_reduced)) {
  chr <- seg_table_reduced[i,]$chr
  family <- seg_table_reduced[i,]$family
  index <- intersect(which(visual.data2$chr == chr), which(visual.data2$family == family))
  index_out <- intersect( index, which(visual.data2$relationship %in% c("c1", "c2", "C")) ) 
  index_in <- intersect( index, which(visual.data2$relationship %in% c("A", "B")) ) 
  rows_to_retrieve <- append(rows_to_retrieve, index_in)
  rows_to_retrieve_compliment <- append(rows_to_retrieve_compliment, index_out)
}


#get rows of data related to MN events
visual.data.results <- cbind(visual.data2[rows_to_retrieve,], infer_MN = rep("yesMN", nrow(visual.data2[rows_to_retrieve,])))
visual.data.results.compliment <- cbind(visual.data2[rows_to_retrieve_compliment,], infer_MN = rep("noMN", nrow(visual.data2[rows_to_retrieve_compliment,])))
visual.data.results.all <- rbind(visual.data.results, visual.data.results.compliment)

#Plot the boxplots and strip plots
boxplot.visual <- ggplot() + 
  #change theme
  theme_bw() + 
  theme(aspect.ratio = 1/4, text=element_text(size=14), axis.text = element_text(color="black"), axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #MN experimental samples 
  #geom_point(data = visual.data.results.all, mapping = aes(x = state, y = TPM, alpha = .2), width = 0.1, position=position_jitterdodge()) +
  geom_boxplot(data = visual.data.results.all, mapping = aes(x = state, y = TPM, alpha = .5), outlier.alpha = 0)  +
  #geom_jitter(data = visual.data.results.all, mapping = aes(x = state, y = TPM, color = infer_MN, alpha = .2), width = 0.1)  +
  #selected samples
  geom_boxplot(data = boxplots.vals, mapping = aes(x = Group, y = TPM, alpha = .5), outlier.alpha = 0)  +
  #all samples
  #geom_jitter(data = visual.data.results.all, mapping = aes(x = "all", y = TPM, alpha = 1), width = 0.1) +
  ylim(c(.25,1.75)) + scale_x_discrete(limits = c("low", "loss", "monosomy", "intermediate_12", "control",  "disomy", "intermediate_32", "gain", "trisomy", "high"), labels=c("low", "ref monosomy", "monosomy", "intermediate_12", "control",  "disomy", "intermediate_32", "ref trisomy", "trisomy", "high")) +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Distribution"))

ggsave(plot = boxplot.visual, filename = "Global_Classifications.pdf", path = sprintf("%s/visual_results", dirpath), device = "pdf", width = 15, height = 4, dpi = 300, units = "in")
