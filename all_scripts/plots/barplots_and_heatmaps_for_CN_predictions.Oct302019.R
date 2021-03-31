#barplots_and_heatmaps_for_CN_predictions.Oct302019.R

#args <- c("/pellmanlab/nikos/Stam_Etai_Scripts", "/pellmanlab/stam_niko/data/processed_bam", "/pellmanlab/nikos/Stam_Etai_Data")
args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
dirpath <- args[2]
datadir <- args[3]
mk_vis_dir <- sprintf("mkdir -p %s/visual_results", dirpath)
system(mk_vis_dir)

require(readxl)
source(sprintf("%s/plots/plot_barplots_of_AllelicAndExpBiasPerSamples.R", scriptsdir))

###########################
# Input variable and data #
###########################

#Cell Annotation list
anno <- data.table(readRDS(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_210129.rds", datadir) ))

#SSM, MA, and Allelic CN information
Ms.chr <- readRDS(file = sprintf("%s/CN_data/CN_predictions.bychr.rds", dirpath))

#Total expression:
chr.TE <- copy(Ms.chr$SSM.TE)
setnames(chr.TE, old = "seqnames", new = "bin_id")
chr.TE <- cbind(chr.TE[,c("bin_id")], 2^(chr.TE[,-c("bin_id")]))
chr.TE <- as.data.frame(chr.TE)

#Allele expression (option B (raw aggs by gene)):
chr.Af <- copy(Ms.chr$raw.A)
setnames(chr.Af, old = "seqnames", new = "bin_id")
chr.Bf <- copy(Ms.chr$raw.B)
setnames(chr.Bf, old = "seqnames", new = "bin_id")
chr.Af <- cbind(chr.Af[,1], chr.Af[,-1]/(chr.Af[, -1] + chr.Bf[, -1]))
chr.Bf <- cbind(chr.Bf[,1], chr.Bf[,-1]/(chr.Af[, -1] + chr.Bf[, -1]))
IDs <- intersect(colnames(chr.TE)[-1], colnames(chr.Af[-1]))

#by sample
plotsdir <- sprintf("%s/visual_results", dirpath)

make_path <- sprintf("mkdir -p %s/bysample/", plotsdir)
system(make_path)

require(doParallel)
require(foreach)
registerDoParallel(26)

foreach(myid = IDs) %dopar% {
  message(myid)
  message(sprintf("Working on %s/bysample/%s.barplots.pdf", plotsdir, myid))
  pdf(file = sprintf("%s/bysample/%s.barplots.pdf", plotsdir, myid), width = 16, height = 12)
  p <- plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = chr.TE, 
                                                     wt = rep(1:nrow(chr.TE)), fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                                     ids = c(myid), chr = "")
  #print(p)
  
  dev.off()
}

#by family:
setkey(anno, WTA.plate)
families <- sort(table(anno[IDs][Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)

make_path <- sprintf("mkdir -p %s/byfamily/", plotsdir)
system(make_path)

require(doParallel)
registerDoParallel(26)

foreach(myfamily = names(families)) %dopar% {
  message(myfamily)
  myids <- anno[Pairs %in% myfamily]$WTA.plate
  pdf(file = sprintf("%s/byfamily/%s.barplots.pdf", plotsdir, myfamily), width = 18, height = 12)
  p <- plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = chr.TE, 
                                                    wt = rep(1:nrow(chr.TE)), fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                                    ids = myids, chr = "")
  #print(p)
  dev.off()
}

#cd ~/Data/ExperimentalData/Papers/MN/reports/CNV12.2/
#rsync -Pvr ejacob@lynx.dfci.harvard.edu:/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/plots.Oct302019 ./
  
#by sample - for alex
plotsdir <- "/pellmanlab/stam_niko/etai_code/experiments_etai/Final_Visualizations_barplots"

make_path <- sprintf("mkdir -p %s/foralex/", plotsdir)
system(make_path)

require(doParallel)
registerDoParallel(13)
IDs <- anno[control_group == "bunch_of_grapes"]$WTA.plate

foreach(myid = IDs) %dopar% {
  message(myid)
  pdf(file = sprintf("%s/foralex/%s.barplots.pdf", plotsdir, myid), width = 16, height = 12)
  p <- plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = chr.TE, 
                                                    wt = rep(1:nrow(chr.TE)), fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                                    ids = c(myid), chr = "")
  print(p)
  
  dev.off()
}

#cd ~/Data/ExperimentalData/Papers/MN/reports/CNV12.2/
#rsync -Pvr ejacob@lynx.dfci.harvard.edu:/singlecellcenter/etai/ExperimentsData/Papers/MN/reports/CNV12.2/Alex_Spector ./
  
################## till here #########################

#allele expression:
#option A (SSM1)
# chr.Af <- copy(Ms.chr$SSM.TEAB.Af)
# setnames(chr.Af, old = "seqnames", new = "bin_id")
# chr.Bf <- copy(Ms.chr$SSM.TEAB.Bf)
# setnames(chr.Bf, old = "seqnames", new = "bin_id")


myids <- c("170512_B8", "170512_B4", "170512_B2")
myids <- triIDs[1:3]
myids <- c("170126_A3", "170126_A4", "170126_A5")
plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = chr.TE, 
                                             wt = rep(1:nrow(chr.TE)), fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                             ids = c(myids), chr = "")




pdf(file = "/singlecellcenter/etai/ExperimentsData/Stamatis/reports.CNV12.0/allPairs.barplots.CNV.v12.0.pdf", width = 15, height = 12)

for(mypair in mypairs) {
  message(mypair)
  myids <- anno[Pairs == mypair]$WTA.plate
  myids <- myids[ myids %in% IDs]
  if(length(myids) < 1) {
    message("No ids..")
    next
  }
  print(myids)
  
  p <- plotAllelicAndExpBiasPerSamples.v8(dfs = chr.TE, 
                                          wt = rep(1:nrow(chr.TE)), fracs = list(Af = chr.Af, Bf = chr.Bf), nogeno = c(),
                                          ids = c(myids), chr = "")
  par(mfrow=c(1,1))
  p <- p + annotation_custom(tableGrob(QC[id %in% myids], rows=NULL), xmin = 20, ymin = max(chr.TE[,myids])*2-0.1,  xmax=Inf, ymax=Inf)
  print(p)
}
dev.off()
