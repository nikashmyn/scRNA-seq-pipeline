#prepare final Results

###############################
### Read in annotation list ###
###############################

anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))

rows1 <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
rows2 <- intersect(which(anno$key_pairs == "Gen1"), which(anno$exclusion == ""))
rows <- union(rows1, rows2)
anno_gen_1_2 <- anno[rows,]

##################################
### Read in possible MN events ###
##################################

#All chromosomes in gen1 and gen2
seg_table <- data.table(read.csv(file = sprintf("%s/aggregated_results/all_chr_profiles.csv", dirpath)))

#All the families are represented in these matrices
possible_events <- data.table(read.csv(file = sprintf("%s/aggregated_results/possible_MN_events.csv", dirpath)))
seg_table_noevent <- data.table(read.csv(file = sprintf("%s/aggregated_results/cells_noMN_events.csv", dirpath)))
(length(unique(seg_table_noevent$family_ids)) + length(unique(possible_events$family_ids))) == length(unique(anno_gen_1_2$Family_IDs))

#all chromosomes are represented in these matrices
seg_table_allevents <- data.table(read.csv(file = sprintf("%s/aggregated_results/all_events.csv", dirpath)))
seg_table_diploid_events <- data.table(read.csv(file = sprintf("%s/aggregated_results/diploid_events.csv", dirpath)))
nrow(seg_table_allevents) + nrow(seg_table_diploid_events) == nrow(seg_table)

#diploid families
seg_table_diploid_families_red <- data.table(read.csv(file = sprintf("%s/aggregated_results/diploid_families.csv", dirpath)))

#all families
all_family_events <- data.table(read.csv(file = sprintf("%s/aggregated_results/all_family_events.csv", dirpath)))

#########################################################
### Create anno lists for included and excluded cells ###
#########################################################

analyzed_cells <- anno[which(anno$Pairs %in% possible_events$family)]
excluded_cells <- anno[which(!anno$Pairs %in% possible_events$family)]
write.csv(analyzed_cells, sprintf("%s/aggregated_results/annotation_list_included.csv", dirpath), row.names=FALSE)
write.csv(excluded_cells, sprintf("%s/aggregated_results/annotation_list_excluded.csv", dirpath), row.names=FALSE)

