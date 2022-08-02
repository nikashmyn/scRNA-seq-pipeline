#visual generator script

##########################
### Script Explanation ###
##########################

#-------------------------------------------------------
#This script was written by Nikos Mynhier 

#####################
### Load Packages ###
#####################

source(sprintf("%s/plots/plot_barplots_of_AllelicAndExpBiasPerSamples4.R", scriptsdir))
source(sprintf("%s/plots/TPM_ratio_func_of_gene_expr.R", scriptsdir))
source(sprintf("%s/plots/plot_TPM_and_ASE_vis.R", scriptsdir))
source(sprintf("%s/plots/plot_AS_TPM_ratio_on_reference.R", scriptsdir))


####################
### loading data ###
####################

##########################Data for Raw Plots###########################
message("Loading data...")

#read in main data
TPM <- readRDS(file=sprintf("%s/aggregated_results/TPM.nolim.rds", dirpath))
TPM_bybin <- readRDS(file=sprintf("%s/aggregated_results/TPM.inv_var.bybin.rds", dirpath))
AS_TPM_bybin <- readRDS(file=sprintf("%s/aggregated_results/AS-TPM.inv_var.bybin.rds", dirpath))
TPM_bychr <- readRDS(file=sprintf("%s/aggregated_results/TPM.inv_var.bychr.rds", dirpath))
ASE_bychr <- readRDS(file=sprintf("%s/aggregated_results/ASE.inv_var.bychr.rds", dirpath))
AS_TPM_bychr <- readRDS(file=sprintf("%s/aggregated_results/AS-TPM.inv_var.bychr.rds", dirpath))

#load cell annotations
anno <- data.table(read.csv( sprintf("%s/work_in_progress/annotation_list.csv", datadir)))
changeCols <- colnames(anno)[-c(13:16)]
anno[,(changeCols):= lapply(.SD, as.character), .SDcols = changeCols] 
anno_red <- anno[which(anno$exclusion == ""),]
columns <- colnames(TPM)[-c(1:6)]
samples_to_visualize <- c(intersect(columns, anno_red$WTA.plate))

#Gene and cell annotaions
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))
centromeres <- read_csv(sprintf("%s/centromeres.csv", datadir))

#read in family events
all_family_events <- read_csv(file = sprintf("%s/aggregated_results/all_family_events.csv", dirpath))

#read in reference distributions
ref_TPM <- readRDS(file = sprintf("%s/aggregated_results/ref_TPM.rds", dirpath))
ref_hap_spec_CN <- readRDS(file = sprintf("%s/aggregated_results/ref_hap_spec_CN.rds", dirpath))

########################################
### Run Individual Plots in Parallel ###
########################################

require(doParallel)
registerDoParallel(10)

#Variables
destDir <- sprintf("%s/visual_results", dirpath)
system(sprintf("mkdir -p %s", destDir)) #make the visual results directory
system(sprintf("mkdir -p %s/individual_plots", destDir)) #make the individual plots directory
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10a", "chr10b", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

#make raw plot pdfs for every sample
foreach(i = c(1:length(samples_to_visualize))) %dopar% {
  myid <- samples_to_visualize[i]
  myfamily_file <- anno[which(anno$WTA.plate == myid),]$Family_IDs
  myid_file <- anno[which(anno$WTA.plate == myid),]$Cell_IDs
  message(myfamily_file)
  #plot barchart for single cell
  plot_barplots_of_AllelicAndExpBiasPerSamples(AS_TPM_bychr = AS_TPM_bychr, myfamily_file = fam_id_name, destDir = sprintf("%s/individual_plots", destDir),
                                               ids = myid, name_ids = myid_file, chr = "", plotAxisText = F, return_plot = T, save_plot = T, single_sample = T)
  #Plot raw plots for relevant chromosomes
  for(chr in chrs) {
    #Chr information
    message(chr)
    #Plot raw data visuals function
    plot_TPM_and_ASE(myfamily_file = myfamily_file, myid_file = myid_file, myid = myid, chr = chr, TPM = TPM_bybin, ASE = AS_TPM_bybin, destDir = sprintf("%s/individual_plots", destDir), controlSampleIDs = controlSampleIDs, combined = T)
    CumulativeTPM_and_ratioTPM_plots(adt = TPM, geneRanges = geneRanges, controlSampleIDs = controlSampleIDs, myfamily_file = myfamily_file, myid_file = myid_file, myid = myid, chr = chr, destDir = sprintf("%s/individual_plots", destDir), binsize = 10)
  }
}


#Make barcharts for each family
setkey(anno, WTA.plate)
families <- sort(table(anno_red[samples_to_visualize][Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }
exclude_chrs <- c()

make_path <- sprintf("mkdir -p %s/byfamily/", destDir)
system(make_path)

foreach(i = c(1:length(names(families)))) %dopar% {
  myfamily = names(families)[i]
  myids <- anno[Pairs %in% myfamily]$WTA.plate
  myids_name <- anno[Pairs %in% myfamily]$Cell_IDs
  fam_id_name <- unique(anno$Family_IDs[which(anno$Pairs == myfamily)])
  message(myfamily) 
  #plot Barcharts by family
  plot_barplots_of_AllelicAndExpBiasPerSamples(AS_TPM_bychr = AS_TPM_bychr, myfamily_file = fam_id_name, destDir = sprintf("%s/byfamily", destDir),
                                                    ids = myids, name_ids = myids_name, chr = "", plotAxisText = F, return_plot = T, save_plot = T)
}

############################################
### Run Family Visual Report in Parallel ###
############################################

#Reduce anno list to relevant gen2 samples
rows1 <- intersect(which(anno$key_pairs == "Gen2"), which(anno$exclusion == ""))
rows2 <- intersect(which(anno$key_pairs == "Gen1"), which(anno$exclusion == ""))
rows <- union(rows1, rows2)
anno_gen_1_2 <- anno[rows,]

#extract family information
setkey(anno_gen_1_2, WTA.plate)
families <- sort(table(anno_gen_1_2[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }

make_path <- sprintf("mkdir -p %s/byfamily/", destDir)
system(make_path)

#Make all cell summary read out by family
#foreach(i = c(1:length(names(families)))) %dopar% { #if not all print its because I start at 10 when last ran
for(i in 1:length(names(families))){
  myfamily = names(families)[i]
  myfamily_file1 = str_remove(myfamily, fixed("(")); myfamily_file = str_remove(myfamily_file1, fixed(")"))
  fam_id_name <- unique(anno$Family_IDs[which(anno$Pairs == myfamily)])
  message(myfamily_file)
  myids <- anno[Pairs %in% myfamily]$WTA.plate
  myids_name <- anno[Pairs %in% myfamily]$Cell_IDs
  pdf(file = sprintf("%s/byfamily/%s.family.pdf", destDir, fam_id_name), width = 18, height = 12)
  #get events for this family
  columns_to_show <- c("generation", "family_ids", "MN_Cell", "MN_Sister", "rupt_time", "chr", "A", "B", "c1", "c2", "CN_pattern", "MN_hap", "MNhap_CN_pattern", "CN_state")
  possible_events_myfam <- all_family_events[which(all_family_events$family == myfamily),columns_to_show]
  #plot table of possible events
  if (nrow(possible_events_myfam) >= 1) {grid.table(possible_events_myfam)}
  #plot Barcharts by family
  q <- plot_barplots_of_AllelicAndExpBiasPerSamples(AS_TPM_bychr = AS_TPM_bychr, myfamily_file = fam_id_name, destDir = sprintf("%s/byfamily", destDir),
                                                    ids = myids, name_ids = myids_name, chr = "", plotAxisText = F, return_plot = T, save_plot = F)
  #Get MN related annotations and do fam specific reference plot 
  MN_related <- anno[union(intersect(which(anno$Pairs == myfamily), which(anno$Sister1 == 1)), intersect(which(anno$Pairs == myfamily), which(anno$Sister2 == 1)))]$WTA.plate
  MN_unrelated <- anno[union(intersect(which(anno$Pairs == myfamily), which(anno$Cousin1 == 1)), intersect(which(anno$Pairs == myfamily), which(anno$Cousin2 == 1)))]$WTA.plate
  plot_AS_TPM_ratio_on_reference(ref_hap_spec_CN = ref_hap_spec_CN, ref_TPM = ref_TPM, AS_TPM_bychr = AS_TPM_bychr, destDir = sprintf("%s/byfamily", destDir), myfamily_file = fam_id_name, myid_files = myids_name, myids = myids, MN_related = MN_related, MN_unrelated = MN_unrelated, return.vis = T, save.vis = T)
  #plot information for possible events
  if (nrow(possible_events_myfam) >= 1 ) {
    if (!is.na(possible_events_myfam$chr)) {
      #Plot raw plots for relevant chromosomes
      for(chr in unique(possible_events_myfam$chr)) {
        #Chr information
        chr <- paste0("chr", chr)
        if(chr == "chr23") {chr <- "chrX"}
        if(chr == "chr10.1") {chr <- "chr10a"}
        if(chr == "chr10.2") {chr <- "chr10b"}
        message(chr)
        files <- data.table()
        files2 <- data.table()
        for(myid in myids[myids %in% samples_to_visualize]) {
          myid_file <- myids_name[which(myids == myid)]
          #Plot raw data visuals function
          files <- rbind(files, plot_TPM_and_ASE(myfamily_file = myfamily_file, myid_file = myid_file, myid = myid, chr = chr, TPM = TPM_bybin, ASE = AS_TPM_bybin, destDir = destDir, controlSampleIDs = controlSampleIDs, save.rds = T, doPlot = F))
          files2 <- rbind(files2, CumulativeTPM_and_ratioTPM_plots(adt = TPM, geneRanges = geneRanges, controlSampleIDs = controlSampleIDs, myfamily_file = myfamily_file, myid_file = myid_file, myid = myid, chr = chr, destDir = destDir, binsize = 10, save.rds = T))
        }
        #TPMratio by bin and var by bin raw plots
        TPM_vis <- lapply(files$TPM_file , function(x){readRDS(x)} )
        var_vis <- lapply(files$var_file , function(x){readRDS(x)} )
        interweave <- c()
        for(i in 1:length(TPM_vis)){interweave <- append(interweave, TPM_vis[i]); interweave <- append(interweave, var_vis[i])}
        do.call("grid.arrange", c(interweave, nrow = length(TPM_vis), ncol = 2)) 
        #CumulativeTPM and ratio_to_expr plots
        expr_vis <- lapply(files2$expr_file , function(x){readRDS(x)} )
        cum_vis <- lapply(files2$cum_file , function(x){readRDS(x)} )
        interweave <- c()
        for(i in 1:length(expr_vis)){interweave <- append(interweave, expr_vis[i]); interweave <- append(interweave, cum_vis[i])}
        do.call("grid.arrange", c(interweave, nrow = length(expr_vis), ncol = 2)) 
      }
    }
  }
  dev.off()
}

#Sanity check that the code ran to completion
print("Done with making visuals")

###########################################################################################
### With plots for each family and all_family_events we begin manual review of the data ###
###########################################################################################
