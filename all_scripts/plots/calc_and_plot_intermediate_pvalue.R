
#function requires the interstat and interstat.arm to be loaded to the workspace prior to calling
calc_and_plot_intermediate_pvalue.arm <- function(myid = "170126_A5", myseqname = "1q", plotme = F, controlSampleIDs) {
  frac <- interstat.arm[sample_id == myid & arm == myseqname]$intermediate_frac
  
  sampfracs <- interstat.arm[sample_id == myid]$intermediate_frac
  sampfracs <- data.table(seqnames = unique(interstat.arm[sample_id == myid]$arm), 
                          fracs = sampfracs,
                          pval = sapply(unique(interstat.arm[sample_id == myid]$arm), 
                                            function(sn) sum(interstat.arm[ sample_id %in% controlSampleIDs & arm == sn]$intermediate_frac > interstat.arm[sample_id == myid  & arm == sn]$intermediate_frac) / 
                                              nrow(interstat.arm[ sample_id %in% controlSampleIDs  & arm == sn])))
  
  
  pval.chr <- sum(interstat.arm[ sample_id %in% controlSampleIDs & arm == myseqname]$intermediate_frac > frac) / 
    nrow(interstat.arm[ sample_id %in% controlSampleIDs  & arm == myseqname])
  
  pval.all <- sum(interstat.arm[ sample_id %in% controlIDs]$intermediate_frac > interstat.arm[sample_id == myid & arm == myseqname]$intermediate_frac) / 
    nrow(interstat.arm[ sample_id %in% controlIDs])
  
  if(plotme) {
    require(ggplot2)
    f <- function(x) {
      r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
      names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
      r
    }
    
    
    p1 <- ggplot(data = interstat.arm[sample_id %in% controlSampleIDs], aes(x = arm, y = intermediate_frac))  + ylab("Fraction of genes") + theme_bw() +
      
      #geom_vline(xintercept = as.numeric(which(unique(interstat.arm$seqnames) == myseqname)), linetype="solid", color = alpha("gray", alpha = 0.4), size = 4) +
      stat_summary(fun.data = f, geom="boxplot", position = position_dodge(2), fill = alpha("lightblue", 0.4), 
                   colour = "black", width = 0.6, size = 0.8, outlier.shape = NA, show.legend = F, alpha=0.6) + #+ ggtitle(sprintf("%s - %s\nMachine Learning Predictions", myid, myseqname)) +
      geom_hline(yintercept=frac, linetype="dotted", color = "red", size = c(2))  + xlab("") + ggtitle(myid) +
      #scale_x_discrete(breaks=posvalues, labels = as.character(round(posvalues/1e6))) + xlab(label = sprintf("Position (Mbp) %s", myseqname)) +
      theme(text = element_text(size=15), axis.text.x = element_text(angle=90, hjust=1)) +
      geom_point(show.legend = F, stat='identity', size=6, 
                 data = sampfracs[seqnames == myseqname],
                 mapping = aes(x=seqnames , y=fracs, fill = "green" ), colour = "black", pch=21)
    
    # do it
    #p1 <- ggplot(data = interstat.arm[sample_id %in% controlSampleIDs], aes(x = seqnames, y = intermediate_frac)) + stat_summary(fun.data = f, geom="boxplot")
    print(p1)
    # boxplot(intermediate_frac ~ seqnames, data = interstat.arm[sample_id %in% controlSampleIDs], outline=FALSE,
    #         las = 2, col = "lightgreen", lwd=2, cex = 1.3, cex.axis = 1.2, cex.lab = 1.3,
    #         ylab = "Fraction of genes", xlab = "")
    # abline(h = frac, col = "red", lwd = 2.5, lty = "dotted")
    
  }
  return(sampfracs)
  return(data.frame(frac = frac, pval.chr = pval.chr, pval.all = pval.all))
}



calc_and_plot_intermediate_pvalue.chr <- function(myid = "170126_A5", myseqname = "chr1", plotme = F, controlSampleIDs) {
  frac <- interstat.chr[sample_id == myid & seqnames == myseqname]$intermediate_frac
  
  sampfracs <- interstat.chr[sample_id == myid]$intermediate_frac
  sampfracs <- data.table(seqnames = unique(interstat.chr[sample_id == myid]$seqnames), 
                          fracs = sampfracs,
                          pval = sapply(unique(interstat.chr[sample_id == myid]$seqnames), 
                                 function(sn) sum(interstat.chr[ sample_id %in% controlSampleIDs & seqnames == sn]$intermediate_frac > interstat.chr[sample_id == myid  & seqnames == sn]$intermediate_frac) / 
                                      nrow(interstat.chr[ sample_id %in% controlSampleIDs  & seqnames == sn])))
  
  #return(data.table(seqnames = unique(interstat.chr[sample_id == myid]$seqnames), fracs = sampfracs))
  pval.chr <- sum(interstat.chr[ sample_id %in% controlSampleIDs & seqnames == myseqname]$intermediate_frac > frac) / 
    nrow(interstat.chr[ sample_id %in% controlSampleIDs  & seqnames == myseqname])
  
  pval.all <- sum(interstat.chr[ sample_id %in% controlIDs]$intermediate_frac > interstat.chr[sample_id == myid & seqnames == myseqname]$intermediate_frac) / 
    nrow(interstat.chr[ sample_id %in% controlIDs])
  
  if(plotme) {
    require(ggplot2)
    f <- function(x) {
      r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
      names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
      r
    }
    
    
    p1 <- ggplot(data = interstat.chr[sample_id %in% controlSampleIDs], aes(x = seqnames, y = intermediate_frac))  + ylab("Fraction of genes") + theme_bw() +
      
      #geom_vline(xintercept = as.numeric(which(unique(interstat.chr$seqnames) == myseqname)), linetype="solid", color = alpha("gray", alpha = 0.4), size = 4) +
      stat_summary(fun.data = f, geom="boxplot", position = position_dodge(2), fill = alpha("lightblue", 0.4), 
                   colour = "black", width = 0.6, size = 0.8, outlier.shape = NA, show.legend = F, alpha=0.6) + #+ ggtitle(sprintf("%s - %s\nMachine Learning Predictions", myid, myseqname)) +
      geom_hline(yintercept=frac, linetype="dotted", color = "red", size = c(2))  + xlab("") + ggtitle(myid) +
      #scale_x_discrete(breaks=posvalues, labels = as.character(round(posvalues/1e6))) + xlab(label = sprintf("Position (Mbp) %s", myseqname)) +
      theme(text = element_text(size=15), axis.text.x = element_text(angle=90, hjust=1)) +
      geom_point(show.legend = F, stat='identity', size=6, 
                 data = sampfracs[seqnames == myseqname],
                 mapping = aes(x=seqnames , y=fracs, fill = "red" ), colour = "black", pch=21)
    
    # do it
    #p1 <- ggplot(data = interstat.chr[sample_id %in% controlSampleIDs], aes(x = seqnames, y = intermediate_frac)) + stat_summary(fun.data = f, geom="boxplot")
    print(p1)
    # boxplot(intermediate_frac ~ seqnames, data = interstat.chr[sample_id %in% controlSampleIDs], outline=FALSE,
    #         las = 2, col = "lightgreen", lwd=2, cex = 1.3, cex.axis = 1.2, cex.lab = 1.3,
    #         ylab = "Fraction of genes", xlab = "")
    # abline(h = frac, col = "red", lwd = 2.5, lty = "dotted")
    
  }
  return(sampfracs)
  return(data.frame(frac = frac, pval.chr = pval.chr, pval.all = pval.all))
}

