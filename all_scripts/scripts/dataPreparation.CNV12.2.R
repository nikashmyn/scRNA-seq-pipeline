
require(data.table)
source('~/WORK/Papers/MLpaper/R/utilities/convert_preds_table_to_matrix_format.R')

lgbmpreds5 <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#preds at gene level

m <- convert_preds_table_to_matrix_format(preds = lgbmpreds5)
saveRDS(m, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/LGBM.classification.rds")

lgbmpreds6 <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/lgbmpreds.bst.basic.vars6.SSM_0.10b.regression.rds")
m <- convert_preds_table_to_matrix_format(preds = lgbmpreds6, model_type = "regression")
saveRDS(m, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/LGBM.regression.rds")


#preds at arm level:
lgbmp <- lgbmpreds5
lgbmp$pred <- max.col(lgbmpreds5[,7:9])
tmp <- lgbmp[, lapply(.SD, function(x) mean(x)), .SDcols = "pred", by = c("sample_id", "arm")]
armpreds <- dcast(tmp, arm ~ sample_id, value.var = "pred" )
armpreds <- armpreds[mixedorder(arm)]
saveRDS(armpreds, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/armpreds.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")

#preds at chr level:
lgbmp <- lgbmpreds5
lgbmp$pred <- max.col(lgbmpreds5[,7:9])
tmp <- lgbmp[, lapply(.SD, function(x) mean(x)), .SDcols = "pred", by = c("sample_id", "seqnames")]
chrpreds <- dcast(tmp, seqnames ~ sample_id, value.var = "pred" )
chrpreds <- chrpreds[mixedorder(seqnames)]
saveRDS(chrpreds, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/chrpreds.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")


#loading data
col_anno <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/col_anno.Apr122018.rds")
controlIDs <- col_anno[experimental_label %in% c("control_LCM_notreatment", "control_gen2_notreatment", "noco_gen2_noMN_lcm", "noco_gen1_noMN", "p53_gen2_noMN")]$WTA.plate



#pval calculations:
#------------------
lgbmpreds5 <- readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/tmp/lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")

#old way:
#tmp <- apply(lgbmpreds5[,7:9], 1, function(x, th = 0.2) (abs(x[2]-x[1]) < th & x[2] > 0.2 & x[1] > 0.2) | (abs(x[3]-x[2]) < th & x[2] > 0.2 & x[3] > 0.2))
maxprob <- apply(lgbmpreds5[, 7:9], 1, max)
lgbmp <- lgbmpreds5
lgbmp$maxprob <- maxprob
lgbmp$intermediate <- maxprob < 0.67
#chr level:
interstat <- lgbmp[, lapply(.SD, function(x) mean(x)), .SDcols = "intermediate", by = c("sample_id", "seqnames")]
interstat$seqnames <- factor(interstat$seqnames, mixedsort(unique(interstat$seqnames)))
setnames(interstat, old = "intermediate", new = "intermediate_frac")
saveRDS(interstat, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.chrlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
#arm level
interstat.arm <- lgbmp[, lapply(.SD, function(x) mean(x)), .SDcols = "intermediate", by = c("sample_id", "arm")]
interstat.arm$arm <- factor(interstat.arm$arm, mixedsort(unique(interstat.arm$arm)))
setnames(interstat.arm, old = "intermediate", new = "intermediate_frac")
saveRDS(interstat.arm, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.armlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")


#entropy pvals calculations:


#lapply(controlSampleIDs, function(myid) {
  
tmp <- lapply(controlSampleIDs, function(myid) {
                  message(myid)
                  rbindlist(lapply(unique(lgbmpreds5$seqnames), function(chr, myid) {
                              data.table(sample_id = myid, 
                                         seqnames = chr, 
                                         t(as.vector(calculate_fraction_of_cns_for_chr(lgbmpreds5[sample_id == myid & seqnames == chr])))) }, myid)) 
                  })

tmp2 <- rbindlist(tmp)
setnames(tmp2, old = names(tmp2)[-c(1:2)], new = c("1","1.5", "2", "2.5", "3"))
tmp2$entropy <- apply(tmp2[, -c(1:2)], 1, function(x) entropy(as.numeric(x)))
fracstat <- tmp2
saveRDS(fracstat, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/fracstat.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")



#------
# tmp
#------
require(gtools)

# interstat <- rbindlist(lapply(unique(lgbmp$seqnames), function(chr) { 
#   rbindlist(lapply(unique(lgbmp$sample_id), 
#                    function(i) data.table(sample_id = i, seqnames = chr, 
#                                           intermediate_frac = mean(lgbmp[ seqnames == chr & sample_id == i ]$intermediate)))) }))
# 
# interstat$seqnames <- factor(interstat$seqnames, mixedsort(unique(interstat$seqnames)))
# saveRDS(interstat, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.chrlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")
# 
# #arm level:
# interstat.arm <- rbindlist(lapply(unique(lgbmp$arm), function(myarm) { 
#   rbindlist(lapply(high_qc_ids, 
#                    function(i) data.table(sample_id = i, arm = myarm, 
#                                           intermediate_frac = mean(lgbmp[ arm == myarm & sample_id == i ]$intermediate)))) }))
# 
# interstat.arm$arm <- factor(interstat.arm$arm, mixedsort(unique(interstat.arm$arm)))
# saveRDS(interstat.arm, file = "/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/data/interstat.armlevel.lgbmpreds.bst.basic.vars6.SSM_0.10b.classification.rds")







#calc_and_plot_intermediate_pvalue.chr(plotme = T, myid = "170126_A5", controlSampleIDs = controlSampleIDs, myseqname = "chr9")
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "170126_A5" & seqnames == "chr17"])
tmp <- raw_data_and_prediction_boxplots(myid = "181013A_6C", chr = "chr1", adt = adt, nonzeros.zs = nonzeros.zs, coding = coding, MLREG = MLREG)

tmp <- raw_data_and_prediction_boxplots(myid = "181013A_6B", chr = "chr1", adt = adt, nonzeros.zs = nonzeros.zs, coding = coding, MLREG = MLREG)


calc_and_plot_intermediate_pvalue.chr(plotme = T, myid = "181013A_1D", controlSampleIDs = controlSampleIDs, myseqname = "chr1")
calc_and_plot_intermediate_pvalue.chr(plotme = T, myid = "181013A_1E", controlSampleIDs = controlSampleIDs, myseqname = "chr1")
calc_and_plot_intermediate_pvalue.chr(plotme = T, myid = "181013A_1A", controlSampleIDs = controlSampleIDs, myseqname = "chr1")
calc_and_plot_intermediate_pvalue.chr(plotme = T, myid = "181013A_1C", controlSampleIDs = controlSampleIDs, myseqname = "chr1")

intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_1D" & seqnames == "chr1"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_1E" & seqnames == "chr1"])
par(mfrow=c(2,1))
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_1A" & seqnames == "chr2"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_1C" & seqnames == "chr2"])

###############
#tmp

calc_and_plot_intermediate_pvalue.chr(plotme = T, myid = "170126_A3", controlSampleIDs = controlIDs, myseqname = "chr1")
calc_and_plot_intermediate_pvalue.arm(plotme = T, myid = "170126_A3", controlSampleIDs = controlIDs, myseqname = "1q")
calc_and_plot_intermediate_pvalue(plotme = T, myid = "170126_A4", controlSampleIDs = controlIDs)
calc_and_plot_intermediate_pvalue(plotme = T, myid = "170126_A5", controlSampleIDs = controlIDs)

calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_6C", controlSampleIDs = controlIDs)
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_6B", controlSampleIDs = controlIDs)
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_6A", controlSampleIDs = controlIDs)
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_5H", controlSampleIDs = controlIDs)

calc_and_plot_intermediate_pvalue.arm(plotme = T, myid = "181013A_6A", controlSampleIDs = controlSampleIDs, myseqname = "1p")
calc_and_plot_intermediate_pvalue.arm(plotme = T, myid = "181013A_5H", controlSampleIDs = controlSampleIDs, myseqname = "1p")
calc_and_plot_intermediate_pvalue.arm(plotme = T, myid = "181013A_6B", controlSampleIDs = controlSampleIDs, myseqname = "1p")
calc_and_plot_intermediate_pvalue.arm(plotme = T, myid = "181013A_6C", controlSampleIDs = controlSampleIDs, myseqname = "1p")

par(mfrow=c(2,1))
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_6B" & seqnames == "chr1"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_6C" & seqnames == "chr1"])

intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_6C" & seqnames == "chr1"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_6B" & seqnames == "chr7"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_6A" & seqnames == "chr1"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_5H" & seqnames == "chr7"])

calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_1A", controlSampleIDs = controlSampleIDs)
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_1C", controlSampleIDs = controlSampleIDs)
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_1D", controlSampleIDs = controlSampleIDs)
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_1E", controlSampleIDs = controlSampleIDs)

par(mfrow=c(4,1))
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_1A" & seqnames == "chr1"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_1C" & seqnames == "chr1"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_1D" & seqnames == "chr1"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "181013A_1E" & seqnames == "chr1"])



calc_and_plot_intermediate_pvalue(plotme = T, myid = "170223_A5a", controlSampleIDs = controlIDs)
calc_and_plot_intermediate_pvalue(plotme = T, myid = "170223_A6a", controlSampleIDs = controlIDs)

calc_and_plot_intermediate_pvalue(plotme = T, myid = "170223_A5a", controlSampleIDs = controlIDs, chr = "chr8")
calc_and_plot_intermediate_pvalue(plotme = T, myid = "170223_A6a", controlSampleIDs = controlIDs, chr = "chr8")
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "170223_A5a" & seqnames == "chr8"])
intermediate <- plot_cn_prediction_probabilities(lgbmpreds5[sample_id == "170223_A6a" & seqnames == "chr8"])


calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_6C", chr = "chr1")
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_6B", chr = "chr1")
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_6A", chr = "chr1")
calc_and_plot_intermediate_pvalue(plotme = T, myid = "181013A_5H", chr = "chr1")

calc_and_plot_intermediate_pvalue(plotme = T, myid = "180815_2A", chr = "chr8")

