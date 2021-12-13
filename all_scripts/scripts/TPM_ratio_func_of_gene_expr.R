CumulativeTPM_and_ratioTPM_plots <- function(destDir = "/pellmanlab/stam_niko/rerun_6_9_2021/data/visual_results",
                                             myfamily_file = "210503_9A9B9C9D", myid = "210503_9A", chr = "chr5", binsize = 10) {

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
  adt <- dt[seqnames != "chrM" & seqnames != "chrY"]
  
  cols <- c("id","seqnames","start","end","width","strand",myid)
  
  adt_all <- adt[,..cols]
  
  adt_control <- data.table(cbind(adt_all, control_TPM = rowMeans(adt[,..controlSampleIDs])))
  
  setkey(adt_control, "control_TPM")
  
  adt_control <- adt_control[which(adt_control$control_TPM > 0),]
  
  adt_control$ratio <- (adt_control[,..myid] / adt_control[,c("control_TPM")])
  
  adt_control_chr <- adt_control[which(adt_control$seqnames == chr),]
  
  #########################################################
  ### Plot TPM ratio as function of avg gene expression ###
  #########################################################
  
  #bin genes
  num_to_rep <- ceiling(nrow(adt_control_chr)/binsize)
  bin <- rep(c(1:num_to_rep), each=binsize)[1:nrow(adt_control_chr)]
  adt_control_bin <- cbind(adt_control_chr, bin)
  
  adt_bin <- adt_control_bin[,-c(1:6)] %>%
    group_by(bin) %>%
    summarise_all(mean, na.rm = TRUE)
  
  adt_bin$control_TPM_logged <- log2(adt_bin$control_TPM)
  
  p1 <- ggplot(adt_bin, aes(x = control_TPM_logged, y = ratio)) + 
    geom_hline(yintercept=c(0, 0.5, 1.0, 1.5, 2), linetype="dashed", color = "black", size = c(0.7, 0.7, 1.5, 0.7, 0.7)) +
    geom_point(color = "red", aes(x = control_TPM_logged, y = ratio)) + 
    ylim(c(0, 3)) + 
    theme_tufte() +
    coord_cartesian(clip="off") +
    geom_rangeframe(y=seq(0, 3, along.with = adt_bin$ratio)) + 
    theme(aspect.ratio = 1/4, text=element_text(size=18), axis.text = element_text(color="black"), axis.title = element_blank()) 
  
  
  ###############################################
  ### Plot cumulative TPM for sample vs ctrl ###
  ###############################################
  
  adt_control_chr$samp_cumsum <- cumsum(adt_control_chr[,..myid])
  adt_control_chr$ctrl_cumsum <- cumsum(adt_control_chr$control_TPM)
  
  p2 <- ggplot(adt_control_chr, aes(x = ctrl_cumsum, y = samp_cumsum)) + 
    geom_point(color = "red", aes(x = ctrl_cumsum, y = samp_cumsum)) + 
    geom_abline(slope = 1.5, intercept = 0) +
    geom_abline(slope = 1.0, intercept = 0) +
    geom_abline(slope = 0.5, intercept = 0) +
    theme_tufte() +
    coord_cartesian(clip="off") +
    geom_rangeframe() + 
    theme(aspect.ratio = 1, text=element_text(size=18), axis.text = element_text(color="black"), axis.title = element_blank()) 

  #save raw visuals in r format for grid arrangement later
  saveRDS(p1, file = sprintf("%s/%s.%s.%s.TPMratiobyExpr.rds", destDir, myfamily_file, myid, chr))
  saveRDS(p2, file = sprintf("%s/%s.%s.%s.CumulativeTPM.rds", destDir, myfamily_file, myid, chr))
  #save raw visuals as pdf for print use
  ggsave(plot = p1, filename = sprintf("%s/%s.%s.%s.TPMratiobyExpr.rds", destDir, myfamily_file, myid, chr), path = destDir, device = "pdf")
  ggsave(plot = p2, filename = sprintf("%s/%s.%s.%s.CumulativeTPM.rds", destDir, myfamily_file, myid, chr), path = destDir, device = "pdf")

  return( data.table(expr_file = sprintf("%s/%s.%s.%s.TPMratiobyExpr.rds", destDir, myfamily_file, myid, chr), cum_file = sprintf("%s/%s.%s.%s.CumulativeTPM.rds", destDir, myfamily_file, myid, chr)) )


}