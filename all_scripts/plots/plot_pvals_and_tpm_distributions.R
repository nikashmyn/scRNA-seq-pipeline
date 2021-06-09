
plot_pvals_and_tpm_distributions <- function(myid=myid, chr=chr, destDir=destDir) {
  
  ### What data is needed ###
  #This data is made by sourcing "calculating_pvals_byarm_nikos.R" before running this script
  #pval_matrix_loss_byarm, pval_matrix_normal_byarm, pval_matrix_gain_byarm, pval_matrix_control_byarm
  #adt_byarm, OLR_preds_byarm, golden_tpms
  pval_matrix_control_byarm <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_control_byarm.rds", dirpath))
  pval_matrix_gain_byarm <- readRDS(file = sprintf("%s/aggregated_results/pval_matrix_gain_byarm.rds", dirpath))
  
  
  
  #save parameters from previous plot in order to revert after this plot
  par.before <- par()
  
  #add chr column to golden samples by adjusting arm strings
  golden_samples$loss$chr <- gsub('.{1}$', '', golden_samples$loss$arm)
  golden_samples$loss$chr <- paste0("chr", golden_samples$loss$chr)
  golden_samples$normal$chr <- gsub('.{1}$', '', golden_samples$normal$arm)
  golden_samples$normal$chr <- paste0("chr", golden_samples$normal$chr)
  golden_samples$gain$chr <- gsub('.{1}$', '', golden_samples$gain$arm)
  golden_samples$gain$chr <- paste0("chr", golden_samples$gain$chr)
  
  #change to i and j to reduce clutter
  arm <- gsub("chr", "", chr)
  p_arm <- gsub('.{0}$', 'p', arm)
  q_arm <- gsub('.{0}$', 'q', arm)
  arms <- c(p_arm, q_arm)

  ###############
  ### VISUALS ###
  ###############
  
  #main data to be plotted byarm
  boxplots.vals <- c()
  boxplots.vals$loss <- golden_tpms$loss_byarm$tpm; boxplots.vals$normal <- golden_tpms$normal_byarm$tpm; boxplots.vals$gain <- golden_tpms$gain_byarm$tpm;
  boxplots.vals$control <- golden_tpms$control_byarm$tpm; boxplots.vals$control_arm_specific <- golden_tpms$control_byarm[which(golden_tpms$control_byarm$arm %in% arms),]$tpm
  dist_colors <- c("deepskyblue4", "darkorchid4", "darkred", "darkgreen", "lightgreen")
  point_colors <- c("orange", "yellow")
  
  layout(matrix(c(1,2), ncol=2))
  
  for (k in 1:2) {
    
    i = myid 
    j = which(adt_byarm$arm %in% arms[k])
    boxplots.vals$control_arm_specific <- golden_tpms$control_byarm[which(golden_tpms$control_byarm$arm %in% arms[k]),]$tpm
    
    #plot dist of arm wide average TPMs for each CN state as labeled in golden set 
    boxplot( boxplots.vals, xlab="CN State", ylab="Normalized Ratio of Arm Average TPM",
             pch = 21, bg = dist_colors, axes = F, frame.plot = FALSE, col = dist_colors,
             ylim = c(.25, 2.25), xlim = c(.5, 5.5), outline = F)
    #Add a grid
    grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
    #Axes labels
    mtext(side = 1, text = names(boxplots.vals), at = c(1,2,3,4,5),
          col = "grey20", line = 1, cex = 0.9)
    mtext(side = 2, text = c(0.5, 1.0, 1.5, 2.0), at = c(0.5, 1.0, 1.5, 2.0),
          col = "grey20", line = 1, cex = 0.9)
    #A line for the OLR prediction CN ratio and average arm tpm 
    abline(a = as.numeric(as.character(OLR_preds_byarm[j,i]))/2, b=0, col=point_colors[1])
    abline(a = as.numeric(as.character(adt_byarm[j,i])), b=0, col=point_colors[2])
    #Build title
    title_txt <- sprintf("%s | %s", colnames(pval_matrix_normal_byarm[j,..i]), pval_matrix_normal_byarm[j,1])
    title(title_txt)
    #Add p-vals to each dist
    text(x=1, y = 2.0, labels = sprintf("P-Val = %s", signif(pval_matrix_loss_byarm[j,..i], digits = 2)), col="black", cex=.6 )
    text(x=2, y = 2.0, labels = sprintf("P-Val = %s", signif(pval_matrix_normal_byarm[j,..i], digits = 2)), col="black", cex=.6 )
    text(x=3, y = .3, labels = sprintf("P-Val = %s", signif(pval_matrix_gain_byarm[j,..i], digits = 2)), col="black", cex=.6 )
    text(x=4, y = .3, labels = sprintf("P-Val = %s", signif(pval_matrix_control_byarm[j,..i], digits = 2)), col="black", cex=.6 )
    text(x=5, y = .3, labels = sprintf("P-Val = %s", signif(pval_matrix_control_byarm_specific[j,..i], digits = 2)), col="black", cex=.6 )
    #Add number of samples in each dist
    text(x=1, y = median(boxplots.vals$loss)+.02, labels = sprintf("Distr. size = %s", length(boxplots.vals$loss)), col="black", cex=.4 )
    text(x=2, y = median(boxplots.vals$normal)+.02, labels = sprintf("Distr. size = %s", length(boxplots.vals$normal)), col="black", cex=.4 )
    text(x=3, y = median(boxplots.vals$gain)+.02, labels = sprintf("Distr. size = %s", length(boxplots.vals$gain)), col="black", cex=.4 )
    text(x=4, y = median(boxplots.vals$control)+.02, labels = sprintf("Distr. size = %s", length(boxplots.vals$control)), col="black", cex=.4 )
    text(x=5, y = median(boxplots.vals$control_arm_specific)+.02, labels = sprintf("Distr. size = %s", length(boxplots.vals$control_arm_specific)), col="black", cex=.4 )
    par(xpd=F)
    legend("topright", legend = c("CN ratio from OLR Prediction", title_txt), col = point_colors, lwd = .5,
           horiz = F, cex = .5, pt.cex = .2, bty = "n", inset = c(-.75,-.045)) # "topright"
  }
  par(par.before)
}
