#inputs
args <- c("/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/scripts", "/pellmanlab/stam_niko/data/processed_bam", "SIS1025a",  "SIS1025b",  "SIS1025d", "SIS1025e", "SIS1025f_Lane1" , "SIS1025f_Lane2")
#args <- commandArgs(trailingOnly = TRUE)
print(args)
scriptsdir <- args[1]
wkpath <- args[2]
experiments <- args[3:length(args)]
mk_agg_dir <- sprintf("mkdir -p %s/aggregated_results", wkpath)
system(mk_agg_dir)

#read in packages and source files
require(data.table)
require(readxl)
require(dplyr)
source('/pellmanlab/nikos/Stam_Etai_Scripts/scripts/runRNAseqPipelineFor_BF9VP.utils.R')
source('/pellmanlab/nikos/Stam_Etai_Scripts/scripts/generate_baseline_datasets_for_ML_and_stats_analysis.R')

#read in all annotation tables
anno1 <- data.table(read_excel(path = "/pellmanlab/nikos/Stam_Etai_Data/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx", sheet = "all", na = c("", "NA")))
anno2 <- data.table(read_excel("/singlecellcenter/etai/ExperimentsData/Stamatis/May072018/Stamatis_list_v12_180508.QCRNAv3.withQCRNAv3.xlsx"))
anno3 <- data.table(read_excel("/pellmanlab/nikos/Stam_Etai_Data/Stamatis_list_v14_181103.xlsx"))
anno4 <- data.table(read_excel("/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Stamatis_list_v15_190910_edit.xlsx"))
anno5 <- data.table(readRDS(file = "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.0/data/col_anno.Apr122018.rds"))

#Get names of all columns
col_anno1 <- colnames(anno1)
col_anno2 <- colnames(anno2)
col_anno3 <- colnames(anno3)
col_anno4 <- colnames(anno4)
col_anno5 <- colnames(anno5)

#Get columns that intersect all anno lists
col_anno12 <- intersect(col_anno1, col_anno2)
col_anno123 <- intersect(col_anno12, col_anno3)
col_anno1234 <- intersect(col_anno123, col_anno4)

#get only intersecting columns from all lists
anno1_intersected <- data.frame(anno1)[,col_anno1234]
anno2_intersected <- data.frame(anno2)[,col_anno1234]
anno3_intersected <- data.frame(anno3)[,col_anno1234]
anno4_intersected <- data.frame(anno4)[,col_anno1234]

#Get combine intersected lists
anno12 <- rbind(anno1_intersected, anno2_intersected)
anno123 <- rbind(anno12, anno3_intersected)
anno1234 <- rbind(anno123, anno4_intersected)

#rbind but dont intersect columns
anno12_long <- rbind(anno1, anno2, fill=T)
anno12_long <- anno12_long[!duplicated(anno12_long$WTA.plate, fromLast = T)]
anno123_long <- rbind(anno12_long, anno3, fill=T)
anno1234_long <- rbind(anno123_long, anno4, fill=T)
anno12345_long <- rbind(anno1234_long, anno5, fill=T)

#Remove first index column in long list
anno12345_long <- anno12345_long[,-c(1)]

#Consolidate Data taking the most recent non-NA value for each group
anno12345_long_dense <- anno12345_long %>%
  group_by(WTA.plate) %>%
    summarise_all( funs(last(na.omit(.))) )

#order columns and eliminate dups
anno1234_order <- anno1234[order(anno1234$WTA.plate),]
anno1234_final <- anno1234_order[!duplicated(anno1234_order$WTA.plate, fromLast = T),]

#Correct some errors where name ends with _
anno1234_final$Fastq_files <- sapply( 1:length(anno1234_final$Fastq_files), function (x) { gsub(pattern="_$", x=anno1234_final$Fastq_files[x], replacement='') } ) 

#Correct some errors where name ends with _ in long list
anno12345_long_dense$Fastq_files <- sapply( 1:length(anno12345_long_dense$Fastq_files), function (x) { gsub(pattern="_$", x=anno12345_long_dense$Fastq_files[x], replacement='') } ) 

#Write out files
write.csv(anno1234_final, "/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Annotation_list_vnikos_210129.csv", row.names = FALSE)
write.csv(anno12345_long_dense, "/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Annotation_list_long_vnikos_210129.csv", row.names = FALSE)

