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
rows <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
anno_gen2 <- anno[rows,]

#reduce anno list to just sisters
rows2 <- union(which(anno_gen2$Sister1 == 1), which(anno_gen2$Sister2 == 1))
anno_gen2 <- anno_gen2[rows2,]

#extract family information
setkey(anno_gen2, WTA.plate)
families <- sort(table(anno_gen2[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }
exclude_chrs <- c()

names <- anno_gen2$WTA.plate

#Prepare Visual Objects
visual.data <- data.table(TPM = as.numeric(flatten(adt.bychr[,..names])))
visual.data$chr <- as.numeric(rep(unique(adt.bychr$seqnames), length = length(visual.data$TPM)))

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

#get state classifications
loss_rows <- intersect( intersect( which(visual.data$pval_loss >= .05), which(visual.data$pval_control < .05)), which(visual.data$pval_gain < .05) )
loss_table <- cbind( rep("monosomy", length(loss_rows)), loss_rows)
control_rows <- intersect( intersect( which(visual.data$pval_loss < .05), which(visual.data$pval_control >= .05)), which(visual.data$pval_gain < .05) )
control_table <- cbind( rep("disomy", length(control_rows)), control_rows)
gain_rows <- intersect( intersect( which(visual.data$pval_loss < .05), which(visual.data$pval_control < .05)), which(visual.data$pval_gain >= .05) )
gain_table <- cbind( rep("trisomy", length(gain_rows)), gain_rows)
ploidy_table <- rbind(rbind(loss_table, control_table), gain_table)
ploidy_table <- data.table(ploidy_table[,c(2,1)]); colnames(ploidy_table) <- c("row", "state")
ploidy_table$row <- as.numeric(ploidy_table$row)
setkey(ploidy_table, "row")

#fill intermediate states
index <- data.table(row = (1:length(visual.data$TPM)))
setkey(index, "row")
ploidy_table2 <- merge(index, ploidy_table, all.x = T) 
ploidy_table2$state[is.na(ploidy_table2$state)] <- "intermediate"

#combine state annotations
visual.data2 <- cbind(visual.data, state = as.factor(ploidy_table2$state))

#remove intermediate classifications
visual.data3 <- visual.data2[which(visual.data2$state != "intermediate"),]

#get intermediate classifications and split into above and below control groups
visual.data.intermediate <- visual.data2[which(visual.data2$state == "intermediate"),]
visual.data.intermediate_low <- visual.data.intermediate[intersect(which(visual.data.intermediate$TPM < 1), which(visual.data.intermediate$TPM > .65)),]
visual.data.intermediate_high <- visual.data.intermediate[intersect(which(visual.data.intermediate$TPM > 1), which(visual.data.intermediate$TPM < 1.35)),]
visual.data.low <- visual.data.intermediate[which(visual.data.intermediate$TPM < .65),]
visual.data.high <- visual.data.intermediate[which(visual.data.intermediate$TPM > 1.35),]


#Plot the boxplots and strip plots
boxplot.visual <- ggplot() + 
  #selected samples
  geom_boxplot(data = boxplots.vals, mapping = aes(x = Group, y = TPM, color = Group), outlier.alpha = 0)  +
  #classified experimental samples 
  geom_jitter(data = visual.data3, mapping = aes(x = state, y = TPM, alpha = .05), width = 0.3)  +
  geom_boxplot(data = visual.data3, mapping = aes(x = state, y = TPM, color = state), outlier.alpha = 0)  +
  #low intermediates
  geom_jitter(data = visual.data.intermediate_low, mapping = aes(x = "low intermediate", y = TPM, alpha = .05), width = 0.3) +
  geom_boxplot(data = visual.data.intermediate_low, mapping = aes(x = "low intermediate", y = TPM, color = "low intermediate"), outlier.alpha = 0) +
  #high intermediates
  geom_jitter(data = visual.data.intermediate_high, mapping = aes(x = "high intermediate", y = TPM, alpha = .05), width = 0.3) +
  geom_boxplot(data = visual.data.intermediate_high, mapping = aes(x = "high intermediate", y = TPM, color = "high intermediate"), outlier.alpha = 0) +
  #very low
  geom_jitter(data = visual.data.low, mapping = aes(x = "very low", y = TPM, alpha = .05), width = 0.3) +
  geom_boxplot(data = visual.data.low, mapping = aes(x = "very low", y = TPM, color = "very low"), outlier.alpha = 0) +
  #very high
  #geom_jitter(data = visual.data.high, mapping = aes(x = "very high", y = TPM, alpha = .05), width = 0.3) +
  #geom_boxplot(data = visual.data.high, mapping = aes(x = "very high", y = TPM, color = "very high"), outlier.alpha = 0) +
  #all samples
  geom_jitter(data = visual.data, mapping = aes(x = "all", y = TPM), width = 0.3) +
  ylim(c(.25,1.75)) + scale_x_discrete(limits = c("loss", "control", "gain", "very low", "monosomy", "low intermediate", "disomy", "high intermediate", "trisomy", "very high", "all"), labels=c("loss", "control", "gain", "very low", "monosomy", "low intermediate", "disomy", "high intermediate", "trisomy", "very high", "all")) +
  labs(x = "Group", y = "TPM Ratio", title = sprintf("TPM Ratio Classes | %s ", myfamily))
