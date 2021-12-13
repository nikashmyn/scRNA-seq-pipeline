source(sprintf('%s/scripts/fromRawTPMtoExprsRatio.R', scriptsdir))

seq_list <- seq(5,125,10); cap_list <- c(); num_genes_list <- c()
for (i in seq_list) {
  expr_ratio <- fromRawTPMtoExprsRatio(rsemtpm, geneRanges, controlSampleIDs, max_pos_col = 6, plusOne = 1, 
                                zerosAsNA = F, normBySd = F, doNormalizeByControls = F,
                                maxExp = 1200, quantileOfExp = 0.99, minDetectionLevel = configs$minDetectionLevel, minNumOfSamplesToDetect = i )
  cap <- expr_ratio$quant99
  num_genes <- nrow(expr_ratio$adt)
  cap_list <- append(cap_list, cap)
  num_genes_list <- append(num_genes_list, num_genes)
}
results <- data.table(num_ctrl_cell = seq_list, cap = cap_list, gene_cnt = num_genes_list)
plot(results$num_ctrl_cell, results$gene_cnt)
plot(results$num_ctrl_cell, results$cap)
