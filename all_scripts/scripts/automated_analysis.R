#automated analysis script

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
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", wkpath)
system(mk_agg_dir)

####################
### Read-in data ###
####################

#Read in data
golden_samples <- readRDS(file = sprintf("%s/aggregated_results/golden_set_tpms_bychr.rds", dirpath))

#Read in grouped sample information
hand_samples <- readRDS(sprintf("%s/aggregated_results/grouped_control_aneuploidies.rds", dirpath))

#read in normed TPM by chr
adt.bychr <- readRDS(file=sprintf("%s/aggregated_results/adt.bychr.rds", dirpath))

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
boxplots.vals <- cbind(golden_samples$loss_bychr$tpm, rep(1, length(golden_samples$loss_bychr$tpm))); 
boxplots.vals <- rbind(boxplots.vals, cbind(golden_samples$control_bychr$tpm, rep(2, length(golden_samples$control_bychr$tpm))))
boxplots.vals <- data.table(rbind(boxplots.vals, cbind(golden_samples$gain_bychr$tpm, rep(3, length(golden_samples$gain_bychr$tpm)))))
colnames(boxplots.vals) <- c("TPM", "Group")
boxplots.vals$TPM <- as.numeric(boxplots.vals$TPM)
boxplots.vals$Group <- as.factor(boxplots.vals$Group)
dist_colors <- c(rep(1, length(golden_samples$loss_bychr$tpm)), rep(10, length(golden_samples$control_bychr$tpm)), rep(20, length(golden_samples$gain_bychr$tpm)))
boxplots.vals$color <- c(dist_colors)
     
#######################
### Strip Plot Data ###
#######################
   
#Reduce anno list to relevant gen2 samples
rows <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
anno_gen2 <- anno[rows,]
#plot the chromosomes for one family
#extract family information
setkey(anno_gen2, WTA.plate)
families <- sort(table(anno_gen2[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }
exclude_chrs <- c()
#pick a family and get related cells
myfamily = names(families)[10]
message(myfamily)
myids <- anno_gen2[Pairs %in% myfamily]$WTA.plate
#Prepare Visual Objects
visual.data <- data.table(TPM = as.numeric(flatten(adt.bychr[,..myids])))
visual.data$chr <- as.numeric(rep(unique(adt.bychr$seqnames), length = length(visual.data$TPM)))
visual.data$cell <- rep(myids, each=length(unique(adt.bychr$seqnames))-length(exclude_chrs))

######################
### Visualize Data ###
######################

#Plot the boxplots and strip plots
boxplot.visual <- ggplot() + 
  geom_jitter(data = visual.data, mapping = aes(x = 4, y = TPM, color = chr, shape = cell), width = 0.2) +
  geom_boxplot(data = boxplots.vals, mapping = aes(x = Group, y = TPM, color = color), outlier.alpha = 0)  +
  ylim(c(.25,1.75)) + scale_x_discrete(limits = c("1", "2", "3", "4"), labels=c("loss", "control", "gain", "family")) +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Classes | %s ", myfamily))

###############################
### Automate Categorization ###
###############################

#add pvals to visual data object
visual.data$pval_loss <- as.numeric(flatten(pval_matrix_loss_bychr[,..myids]))
visual.data$pval_control <- as.numeric(flatten(pval_matrix_control_bychr[,..myids]))
visual.data$pval_gain <- as.numeric(flatten(pval_matrix_gain_bychr[,..myids]))

#get data for possible MN events
possible_chromosomes <- unique(visual.data[which(visual.data$pval_control <= .05),]$chr)
possible_data <- visual.data[which(visual.data$chr %in% possible_chromosomes),]

#get and bind sister and cousin information to possible MN event data
Cell_lineage_data <- anno_gen2[which(anno_gen2$WTA.plate %in% unique(possible_data$cell)), c("WTA.plate", "Sister1", "Sister2", "Cousin1", "Cousin2")]
family_info <- cbind("c1", Cell_lineage_data$WTA.plate[which(Cell_lineage_data$Sister1 == 1)])
family_info <- rbind(family_info, c("c2", Cell_lineage_data$WTA.plate[which(Cell_lineage_data$Sister2 == 1)]))
family_info <- rbind(family_info, c("c3", Cell_lineage_data$WTA.plate[which(Cell_lineage_data$Cousin1 == 1)]))
family_info <- rbind(family_info, c("c4", Cell_lineage_data$WTA.plate[which(Cell_lineage_data$Cousin2 == 1)]))
family_info <- data.table(family_info); colnames(family_info) <- c("relationship", "cell")

#merge cell lineage data with possible events
setkey(family_info, "cell"); setkey(possible_data, "cell")
possible_data2 <- merge(possible_data, family_info)
setkey(possible_data2, "chr")
possible_data2 <- possible_data2[,c("chr", "relationship", "cell", "TPM", "pval_loss", "pval_control", "pval_gain")]

#Create state annotation based on pval
loss_rows <- intersect( intersect( which(possible_data2$pval_loss >= .05), which(possible_data2$pval_control < .05)), which(possible_data2$pval_gain < .05) )
loss_table <- cbind( rep("monosomy", length(loss_rows)), loss_rows)
control_rows <- intersect( intersect( which(possible_data2$pval_loss < .05), which(possible_data2$pval_control >= .05)), which(possible_data2$pval_gain < .05) )
control_table <- cbind( rep("disomy", length(control_rows)), control_rows)
gain_rows <- intersect( intersect( which(possible_data2$pval_loss < .05), which(possible_data2$pval_control < .05)), which(possible_data2$pval_gain >= .05) )
gain_table <- cbind( rep("trisomy", length(gain_rows)), gain_rows)
ploidy_table <- rbind(rbind(loss_table, control_table), gain_table)
ploidy_table <- data.table(ploidy_table[,c(2,1)]); colnames(ploidy_table) <- c("row", "state")
ploidy_table$row <- as.numeric(ploidy_table$row)
setkey(ploidy_table, "row")

#fill in holes in the categorization as unknown
index <- data.table(row = (1:length(possible_data2$relationship)))
setkey(index, "row")
ploidy_table2 <- merge(index, ploidy_table, all.x = T) 
ploidy_table2$state[is.na(ploidy_table2$state)] <- "intermediate"

#bind state annotations with possible event data
possible_data3 <- cbind(possible_data2, state = ploidy_table2$state)
keycol <-c("chr","relationship")
setorderv(possible_data3, keycol)

### layout possible scenarios ###
scenarios <- data.table(#full defect 2:1
                        full_defect_21 = c("monosomy", "monosomy", "disomy", "disomy"),
                        #no defect 2:1
                        no_defect_21_a = c("monosomy", "disomy", "disomy", "disomy"),
                        no_defect_21_b = c("disomy", "monosomy", "disomy", "disomy"),
                        #full defect 3:2
                        full_defect_32 = c("disomy", "disomy", "monosomy", "monosomy"),
                        #no defect 3:2
                        no_defect_32_a = c("disomy", "trisomy", "monosomy", "monosomy"),
                        no_defect_32_b = c("trisomy", "disomy", "monosomy", "monosomy"))

#label possible data with these scenarios
possible_data4 <- possible_data3[,c("chr", "relationship", "state")]
possible_data5 <- cbind(possible_data4, scenarios)



