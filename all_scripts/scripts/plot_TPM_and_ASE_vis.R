#plot_TPM_and_ASE_vis.R

plot_TPM_and_ASE <- function(myfamily_file = "F258", myid_file = "F258.3", myid = "210503_9C", chr = "chr5", destDir = "/pellmanlab/stam_niko/rerun_6_9_2021/data/visual_results", TPM, ASE, controlSampleIDs, doPlot = T, save.rds = F) { 
  
  #Get the centromere position
  if(chr %in% c("chr10a", "chr10b")) {chr_tmp <- "chr10"} else {chr_tmp <- chr}
  rng <- range(centromeres[seqnames==chr_tmp, c("start", "end")])
  
  ###################
  ### TPM Section ###
  ###################
  
  #reduce matrix to one chr one cell
  TPM <- data.table(TPM)
  columns <- c("chr", "start", "end", "mid", myid)
  rows <- which(TPM$chr == chr)
  TPM <- TPM[rows,..columns]
  TPM$mid <- (TPM$start + TPM$end) / 2
  setnames(TPM, old = c(myid), new = c("cell"))
  
  #create continuous ticks
  all_ticks <- seq(0,250000000,5000000)
  major_ticks <- seq(0,250000000,25000000)
  all_ticks_breaks <- all_ticks[all_ticks < max(TPM$end)]
  all_ticks_labs <- as.character(round(all_ticks_breaks/1e6))
  all_ticks_labs[which(!all_ticks_breaks %in% major_ticks)] <- ""
  
  #Plot TPM data
  p1 <- ggplot(TPM, aes(x = mid, y = cell)) + 
    geom_vline(xintercept=min(TPM$mid[which(TPM$mid >= rng[1])]), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
    geom_hline(yintercept=c(0, 0.5, 1.0, 1.5, 2), linetype="dashed", color = "black", size = c(0.7, 0.7, 1.5, 0.7, 0.7)) +
    geom_point(color = "black", aes(x = mid, y = cell), size = 3) + 
    scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
    scale_x_continuous(breaks=all_ticks_breaks, labels = all_ticks_labs, limits = c(0, 250000000)) + 
    coord_cartesian(clip="off") +
    geom_rangeframe(x=seq(0, max(all_ticks_breaks), along.with = TPM$mid), y=seq(0, 2.5, along.with = TPM$mid)) + 
    theme_tufte() +
    theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"),
          axis.text = element_text(color="black"), axis.title = element_blank()) 
  
  
  ###################
  ### ASE Section ###
  ###################
  
  #Reduce the matrix to one chromosome in one cell
  ASE$combined <- rbind(cbind(Hap = "A", ASE$A), cbind(Hap = "B", ASE$B))
  columns <- c("Hap", "chr", "start", "end", myid)
  rows <- which(ASE$combined$chr == chr)
  ASE_comb_chr_cell <- ASE$combined[rows,..columns]
  ASE_comb_chr_cell$mid <- (ASE_comb_chr_cell$start + ASE_comb_chr_cell$end) / 2
  setnames(ASE_comb_chr_cell, old = c(myid), new = c("cell"))
  
  #create continuous ticks
  all_ticks <- seq(0,250000000,5000000)
  major_ticks <- seq(0,250000000,25000000)
  all_ticks_breaks <- all_ticks[all_ticks < max(ASE_comb_chr_cell$end)]; 
  all_ticks_labs <- as.character(round(all_ticks_breaks/1e6))
  all_ticks_labs[which(!all_ticks_breaks %in% major_ticks)] <- ""
  
  #Plot the ASE data
  p2 <- ggplot(data = ASE_comb_chr_cell, mapping = aes(x = mid, y = cell, fill=Hap)) + 
    geom_hline(yintercept=c(0, 1.0, 2), linetype="dashed", color = "black", size = c(0.7, 1.5, 0.7)) +
    geom_vline(xintercept=min(ASE_comb_chr_cell$mid[which(ASE_comb_chr_cell$mid >= rng[1])]), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
    geom_point(position = position_dodge(0.7), show.legend = F, stat='identity', size=3, data = ASE_comb_chr_cell, 
               mapping = aes(x=mid, y=cell, color = Hap)) + 
    scale_colour_manual(values = c("dodgerblue3", "firebrick2")) +
    scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
    scale_x_continuous(breaks=all_ticks_breaks, labels = all_ticks_labs, limits = c(0, 250000000)) +
    coord_cartesian(clip="off") +
    geom_rangeframe(x=seq(min(ASE_comb_chr_cell$start), max(ASE_comb_chr_cell$mid), along.with = ASE_comb_chr_cell$mid), y=seq(0, 2.5, along.with = ASE_comb_chr_cell$cell)) + 
    theme_tufte() +
    theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"), 
          axis.text = element_text(color="black"), axis.title = element_blank())
  
  #save the visuals 
  ggsave(plot = p1, filename = sprintf("%s.%s.%s.TotalTPM.pdf", myfamily_file, myid_file, chr), path = destDir, device = "pdf", width = 12, height = 3, dpi = 300, units = "in")
  ggsave(plot = p2, filename = sprintf("%s.%s.%s.AllelicCov.pdf", myfamily_file, myid_file, chr), path = destDir, device = "pdf", width = 12, height = 3, dpi = 300, units = "in")
  if(save.rds == T) {
    saveRDS(p1, file = sprintf("%s/%s.%s.%s.TotalTPM.rds", destDir, myfamily_file, myid_file, chr))
    saveRDS(p2, file = sprintf("%s/%s.%s.%s.AllelicCov.rds", destDir, myfamily_file, myid_file, chr))
    return( data.table(TPM_file = sprintf("%s/%s.%s.%s.TotalTPM.rds", destDir, myfamily_file, myid_file, chr), var_file = sprintf("%s/%s.%s.%s.AllelicCov.rds", destDir, myfamily_file, myid_file, chr)) )
  }
  
  #draw plots if T
  if(doPlot) {
    grid.arrange(p1, p2, nrow = 2)
  }
  
}


