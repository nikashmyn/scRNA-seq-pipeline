############################################################################################
# Input Desired script Directories #
############################################################################################

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/pellmanlab/nikos/scRNA-seq-pipeline/all_scripts", "/pellmanlab/stam_niko/rerun_6_9_2021/data", "/pellmanlab/nikos/Stam_Etai_Data", "SIS1025a",  "SIS1025b", "SIS1025c", "SIS1025d", "SIS1025e", "SIS1025f_Lane1", "SIS1025f_Lane2", "SIS1025g_Lane1", "SIS1025g_Lane2", "SIS1025misc", "SIS1025targ")
print(args)
scriptsdir <- args[1]
wkpath <- dirpath <- args[2]
datadir <- args[3]
experiments <- args[4:length(args)]

############################################################################################
# Read in data #
############################################################################################

#Read in data binned by genomic coordinates
TPM_bybin <- readRDS(file=sprintf("%s/aggregated_results/TPM.inv_var.bybin.rds", dirpath))
ASE_bybin <- readRDS(file=sprintf("%s/aggregated_results/ASE.inv_var.bybin.rds", dirpath))
AS_TPM_bybin <- readRDS(file=sprintf("%s/aggregated_results/AS-TPM.inv_var.bybin.rds", dirpath))

#Control samples
controlSampleIDs <- readRDS(sprintf("%s/aggregated_results/controlSampleIDs.rds", wkpath))

############################################################################################
# Parameters #
############################################################################################

#Pick chromosome to plot
myid = "171205_1C"
chr = "chr10"

#Get the centromere position
centromeres <- readRDS( sprintf("%s/centromeres.rds", datadir) )
rng <- range(centromeres[seqnames==chr, c("start", "end")])

############################################################################################
# ASE Visual #
############################################################################################

##combine A and B data frames
#ASE_bybin$combined <- rbind( cbind(Hap = "A" , ASE_bybin$A) , cbind(Hap = "B" , ASE_bybin$B) )
#
##Reduce the matrix to one chromosome in one cell
#columns <- c("Hap", "chr", "start", "end", myid)
#rows <- which(ASE_bybin$combined$chr == chr)
#ASE_comb_chr_cell <- ASE_bybin$combined[rows,..columns]
#ASE_comb_chr_cell$mid <- (ASE_comb_chr_cell$start + ASE_comb_chr_cell$end) / 2
#setnames(ASE_comb_chr_cell, old = c(myid), new = c("cell"))
#
##create continuous ticks
#all_ticks <- seq(0,250000000,5000000)
#major_ticks <- seq(0,250000000,25000000)
#all_ticks_breaks <- all_ticks[all_ticks < max(ASE_comb_chr_cell$end)]; 
#all_ticks_labs <- as.character(round(all_ticks_breaks/1e6))
#all_ticks_labs[which(!all_ticks_breaks %in% major_ticks)] <- ""
#
#p2 <- ggplot(data = ASE_comb_chr_cell, mapping = aes(x = mid, y = cell, fill=Hap)) + 
#  geom_hline(yintercept=c(0, 1.0, 2), linetype="dashed", color = "black", size = c(0.7, 1.5, 0.7)) +
#  geom_vline(xintercept=min(ASE_comb_chr_cell$mid[which(ASE_comb_chr_cell$mid >= rng[1])]), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
#  geom_point(position = position_dodge(0.7), show.legend = F, stat='identity', size=3, data = ASE_comb_chr_cell, 
#             mapping = aes(x=mid, y=cell, color = Hap)) + 
#  scale_colour_manual(values = c("dodgerblue3", "firebrick2")) +
#  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
#  scale_x_continuous(breaks=all_ticks_breaks, labels = all_ticks_labs, limits = c(0, 250000000)) +
#  coord_cartesian(clip="off") +
#  geom_rangeframe(x=seq(min(ASE_comb_chr_cell$start), max(ASE_comb_chr_cell$mid), along.with = ASE_comb_chr_cell$mid), y=seq(0, 2.5, along.with = ASE_comb_chr_cell$cell)) + 
#  theme_tufte() +
#  theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"), 
#        axis.text = element_text(color="black"), axis.title = element_blank())

############################################################################################
# TPM Section #
############################################################################################

#reduce matrix to one chr one cell
columns <- c("chr", "start", "end", myid)
rows <- which(TPM_bybin$chr == chr)
TPM <- TPM_bybin[rows,..columns]
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

############################################################################################
# AS-TPM Visual #
############################################################################################

#combine A and B data frames
AS_TPM_bybin$combined <- rbind( cbind(Hap = "A" , AS_TPM_bybin$A) , cbind(Hap = "B" , AS_TPM_bybin$B) )

#Reduce the matrix to one chromosome in one cell
columns <- c("Hap", "chr", "start", "end", myid)
rows <- which(AS_TPM_bybin$combined$chr == chr)
AS_TPM_comb_chr_cell <- AS_TPM_bybin$combined[rows,..columns]
AS_TPM_comb_chr_cell$mid <- (AS_TPM_comb_chr_cell$start + AS_TPM_comb_chr_cell$end) / 2
setnames(AS_TPM_comb_chr_cell, old = c(myid), new = c("cell"))

#create continuous ticks
all_ticks <- seq(0,250000000,5000000)
major_ticks <- seq(0,250000000,25000000)
all_ticks_breaks <- all_ticks[all_ticks < max(AS_TPM_comb_chr_cell$end)]; 
all_ticks_labs <- as.character(round(all_ticks_breaks/1e6))
all_ticks_labs[which(!all_ticks_breaks %in% major_ticks)] <- ""

p3 <- ggplot(data = AS_TPM_comb_chr_cell, mapping = aes(x = mid, y = cell, fill=Hap)) + 
  geom_hline(yintercept=c(0, 1.0, 2), linetype="dashed", color = "black", size = c(0.7, 1.5, 0.7)) +
  geom_vline(xintercept=min(AS_TPM_comb_chr_cell$mid[which(AS_TPM_comb_chr_cell$mid >= rng[1])]), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
  geom_point(position = position_dodge(0.7), show.legend = F, stat='identity', size=3, data = AS_TPM_comb_chr_cell, 
             mapping = aes(x=mid, y=cell, color = Hap)) + 
  scale_colour_manual(values = c("dodgerblue3", "firebrick2")) +
  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
  scale_x_continuous(breaks=all_ticks_breaks, labels = all_ticks_labs, limits = c(0, 250000000)) +
  coord_cartesian(clip="off") +
  geom_rangeframe(x=seq(min(AS_TPM_comb_chr_cell$start), max(AS_TPM_comb_chr_cell$mid), along.with = AS_TPM_comb_chr_cell$mid), y=seq(0, 2.5, along.with = AS_TPM_comb_chr_cell$cell)) + 
  theme_tufte() +
  theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"), 
        axis.text = element_text(color="black"), axis.title = element_blank())

#############################################################################################
## AS-TPM v2 Visual #
#############################################################################################
#
##combine A and B data frames
#AS_TPM_normed_v2$combined <- rbind( cbind(Hap = "A" , AS_TPM_normed_v2$A) , cbind(Hap = "B" , AS_TPM_normed_v2$B) )
#
##Reduce the matrix to one chromosome in one cell
#columns <- c("Hap", "chr", "start", "end", myid)
#rows <- which(AS_TPM_normed_v2$combined$chr == chr)
#AS_TPM_comb_chr_cell <- AS_TPM_normed_v2$combined[rows,..columns]
#AS_TPM_comb_chr_cell$mid <- (AS_TPM_comb_chr_cell$start + AS_TPM_comb_chr_cell$end) / 2
#setnames(AS_TPM_comb_chr_cell, old = c(myid), new = c("cell"))
#
##create continuous ticks
#all_ticks <- seq(0,250000000,5000000)
#major_ticks <- seq(0,250000000,25000000)
#all_ticks_breaks <- all_ticks[all_ticks < max(AS_TPM_comb_chr_cell$end)]; 
#all_ticks_labs <- as.character(round(all_ticks_breaks/1e6))
#all_ticks_labs[which(!all_ticks_breaks %in% major_ticks)] <- ""
#
#p4 <- ggplot(data = AS_TPM_comb_chr_cell, mapping = aes(x = mid, y = cell, fill=Hap)) + 
#  geom_hline(yintercept=c(0, 1.0, 2), linetype="dashed", color = "black", size = c(0.7, 1.5, 0.7)) +
#  geom_vline(xintercept=min(AS_TPM_comb_chr_cell$mid[which(AS_TPM_comb_chr_cell$mid >= rng[1])]), linetype="solid", color = alpha("black", alpha = 0.4), size = 5) +
#  geom_point(position = position_dodge(0.7), show.legend = F, stat='identity', size=3, data = AS_TPM_comb_chr_cell, 
#             mapping = aes(x=mid, y=cell, color = Hap)) + 
#  scale_colour_manual(values = c("dodgerblue3", "firebrick2")) +
#  scale_y_continuous(breaks = c(0,.5,1,1.5,2,2.5), labels = c("0","","1","","2",""), limits = c(0,2.5)) +
#  scale_x_continuous(breaks=all_ticks_breaks, labels = all_ticks_labs, limits = c(0, 250000000)) +
#  coord_cartesian(clip="off") +
#  geom_rangeframe(x=seq(min(AS_TPM_comb_chr_cell$start), max(AS_TPM_comb_chr_cell$mid), along.with = AS_TPM_comb_chr_cell$mid), y=seq(0, 2.5, along.with = AS_TPM_comb_chr_cell$cell)) + 
#  theme_tufte() +
#  theme(aspect.ratio = 1/5, text=element_text(size=18), axis.ticks.length=unit(.25, "cm"), 
#        axis.text = element_text(color="black"), axis.title = element_blank())

############################################################################################
# arrange and output combined visual  #
############################################################################################

grid.arrange(p1, p3, ncol = 2)

############################################################################################
# Barchart visual  #
############################################################################################

#read in normed TPM by chr
TPM_bychr <- readRDS(file=sprintf("%s/aggregated_results/TPM.inv_var.bychr.rds", dirpath))
TPM_bychr2 <- cbind( TPM_bychr[,c(1)], (2*TPM_bychr[,-c(1)]) )
setnames(TPM_bychr2, old = "chr", new = "bin_id")


#new ASE
ASE_bychr <- readRDS(file=sprintf("%s/aggregated_results/ASE.inv_var.bychr.rds", dirpath))
VAR_A <- copy(ASE_bychr$AF)
VAR_B <- cbind(ASE_bychr$AF[,c(1)], 1-ASE_bychr$AF[,-c(1)])
setnames(VAR_A, old = "chr", new = "bin_id")
setnames(VAR_B, old = "chr", new = "bin_id")

source(sprintf("%s/plots/plot_barplots_of_AllelicAndExpBiasPerSamples2.R", scriptsdir))

#extract family information
setkey(anno_gen_1_2, WTA.plate)
families <- sort(table(anno_gen_1_2[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)
if ("" %in% names(families)) { families <- families[which(!names(families) == "")] }

destDir <- sprintf("%s/visual_results", dirpath)
system(sprintf("mkdir -p %s", destDir))
make_path <- sprintf("mkdir -p %s/byfamily/", destDir)
system(make_path)

myfamily = names(families)[i]
myids <- anno[Pairs %in% myfamily]$WTA.plate
myids_name <- anno[Pairs %in% myfamily]$Cell_IDs
fam_id_name <- unique(anno$Family_IDs[which(anno$Pairs == myfamily)])
message(myfamily)

q <- plot_barplots_of_AllelicAndExpBiasPerSamples(dfs = TPM_bychr2, myfamily_file = fam_id_name, destDir = sprintf("%s/byfamily", destDir),
                                                  wt = rep(1:nrow(TPM_bychr)), fracs = list(Af = VAR_A, Bf = VAR_B), nogeno = c(),
                                                  ids = myids, name_ids = myids_name, chr = "", plotAxisText = F, return_plot = T)

