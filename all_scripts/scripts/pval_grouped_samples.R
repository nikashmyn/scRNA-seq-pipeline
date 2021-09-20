#########################################
### Pval for manually grouped samples ###
#########################################

#source cluster script
source(sprintf("%s/plots/ClusterMap_byfamily.R", scriptsdir))

#Manual input of grouped samples
hand_sample_ids <- c("180815_1A", "180815_1A", "180815_1B", "180815_1D", "180815_1E", "180815_1F", "180815_1F", "180815_1F", "180815_1G", "180815_1G", "180815_1H", "180815_1H", "180815_2A", "180815_2E", "180815_2F", "180815_2G", "181013A_1D", "181013A_1E", "181013A_2F", "181013A_2F", "181013A_2G", "181013A_2G", "181013A_2H", "181013A_3A", "181013A_5H", "181013A_6A", "190701_9F", "190702_6D", "171205_1D", "171205_3C", "171205_3D", "171205_3D", "190627_4E", "190628_2A", "190628_3C", "190701_3A", "171205_2G", "171205_2H", "180301_1A", "190701_3C", "190701_3D")
hand_sample_CNs <- c(1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 1, 1, 1, 1, 1, 1, 3, 3, 1, 1, 1, 1, 3, 3, 3, 3, 3)
hand_sample_chrs <- c(1, 15, 15, 15, 15, 1, 2, 8, 2, 8, 2, 8, 2, 2, 2, 8, 1, 1, 2, 8, 2, 8, 2, 2, 1, 1, 9, 7, 12, 8, 1, 8, 1, 11, 17, 2, 12, 12, 12, 12, 12)

#hand_sample_ids2 <- c("161130_B6", "161130_B6", "161130_B7", "161130_B7", "161130_B7", "161229_A1", "161229_A1", "161229_A1", "161229_A1", "161229_A1", "161229_A2", "161229_A2", "161229_A2", "161229_A2", "161229_A2", "170119_A5", "170119_A5", "170119_A6", "170119_A6", "170119_A6", "170126_A4", "170126_A5", "170126_A8a", "170126_A8b", "170126_B1a", "180815_3C", "180815_3C", "180815_3C", "180815_3C", "180815_3D", "180815_3D", "180815_3D", "180815_3D", "190702_3C", "190702_3C", "190702_3D", "190702_3D", "190702_3E", "190702_3E")
#hand_sample_CNs2 <- c(1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 1, 3, 3, 3, 3, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 1, 1, 3, 3, 1, 1, 3, 1, 1, 3, 3, 3, 3)
#hand_sample_chrs2 <- c(19, 20, 11, 19, 20, 2, 4, 5, 6, 20, 2, 4, 5, 6, 20, 6, 7, 3, 6, 7, 17, 17, 8, 8, 19, 11, 14, 16, 20, 11, 14, 16, 20, 7, 14, 7, 14, 7, 14 )
#
#hand_sample_ids <- append(hand_sample_ids, hand_sample_ids2)
#hand_sample_CNs <- append(hand_sample_CNs, hand_sample_CNs2)
#hand_sample_chrs <- append(hand_sample_chrs, hand_sample_chrs2)

hand_samples <- data.table(ID = hand_sample_ids, CN = hand_sample_CNs, chr = hand_sample_chrs)
saveRDS(hand_samples, sprintf("%s/aggregated_results/grouped_control_aneuploidies.rds", dirpath))

#myfamily = "Manual Sample Grouping"
#message(myfamily)
#myids <- hand_samples$ID
#mychrs <- hand_samples$chr
#
##Plot cluster map output by family
#ClusterMap3D_Groups(adt.bychr, var_mat_sep.bychr, abs_allele_diff_mat.bychr, myids, myfamily, mychrs)
#
##source pval calculations. v3 uses new manual aneuploid grouping. 
#source(sprintf("%s/scripts/calculating_control_pvals_bychr_v3.R", scriptsdir))
#
##load plotting script and plot pvals from grouped samples
#source(sprintf('%s/plots/plot_pvals_and_tpm_distributions_bychr.R', scriptsdir))
#plot_pvals_and_tpm_distributions(myid=hand_samples$ID[12], chr=paste0("chr", hand_samples$chr[12]) , destDir=destDir)

