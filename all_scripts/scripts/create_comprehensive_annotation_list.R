#inputs
#args <- c("/pellmanlab/stam_niko/etai_code/DFCI.scRNAseq.workflows/scripts", "/pellmanlab/stam_niko/data/processed_bam", "SIS1025a",  "SIS1025b",  "SIS1025d", "SIS1025e", "SIS1025f_Lane1" , "SIS1025f_Lane2")
args <- commandArgs(trailingOnly = TRUE)
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
anno6 <- data.table(read.csv("/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/annotation_list_SIS1025g_nikos.csv"))

#Get names of all columns
col_anno1 <- colnames(anno1)
col_anno2 <- colnames(anno2)
col_anno3 <- colnames(anno3)
col_anno4 <- colnames(anno4)
col_anno5 <- colnames(anno5)
col_anno6 <- colnames(anno6)

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
anno123456_long <- rbind(anno12345_long, anno6, fill=T)

#Remove first index column in long list
anno123456_long <- anno123456_long[,-c(1)]

#Consolidate Data taking the most recent non-NA value for each group
anno123456_long_dense <- anno123456_long %>%
  group_by(WTA.plate) %>%
    summarise_all( funs(last(na.omit(.))) )

#order columns and eliminate dups
anno1234_order <- anno1234[order(anno1234$WTA.plate),]
anno1234_final <- anno1234_order[!duplicated(anno1234_order$WTA.plate, fromLast = T),]

#Correct some errors where name ends with _
anno1234_final$Fastq_files <- sapply( 1:length(anno1234_final$Fastq_files), function (x) { gsub(pattern="_$", x=anno1234_final$Fastq_files[x], replacement='') } ) 

#Correct some errors where name ends with _ in long list
anno123456_long_dense$Fastq_files <- sapply( 1:length(anno123456_long_dense$Fastq_files), function (x) { gsub(pattern="_$", x=anno123456_long_dense$Fastq_files[x], replacement='') } ) 

#Reorder and manually select some columns
columns <- colnames(anno123456_long_dense)[c(1:36,165:168)]
anno123456_long_dense_reduced <- anno123456_long_dense[,columns]

#Write out files
write.csv(anno1234_final, "/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Annotation_list_vnikos_210129.csv", row.names = FALSE)
write.csv(anno123456_long_dense, "/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Annotation_list_long_vnikos_210129.csv", row.names = FALSE)
write.csv(anno123456_long_dense_reduced, "/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Annotation_list_long_vnikos_210129.csv", row.names = FALSE)

#Combine the aggregate list with the corrected annotation list
anno_fixed <- data.table(read.csv( sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir)))
anno123456_long_fixed <- rbind(anno_fixed, anno123456_long_dense, fill=T)
anno123456_long_fixed_dense <- anno123456_long_fixed %>%
  group_by(WTA.plate) %>%
  summarise_all( funs(first(na.omit(.))) )

#Remove extra space in "control " in event column.
anno123456_long_fixed_dense$event[anno123456_long_fixed_dense$event == "control "] <- "control"

#Manually remove unsequenced sample
anno123456_long_fixed_dense <- anno123456_long_fixed_dense[-c(which(anno123456_long_fixed_dense$WTA.plate == "210621_10H")),]

write.csv(anno123456_long_fixed_dense, sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir), row.names = FALSE)



#######################################################################################


#Generate final list by excluding irrelevant samples and columns
anno_all <- data.table(read.csv("/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/SS4_list_Oct13_St.csv", header = F))

colnames(anno_all) <- as.character(anno_all[8,]) #use colnames row as headers
anno_all <- anno_all[,c(1:43)] #get rid of empty columns
anno_all <- anno_all[-c(1:9),] #Get rid of empty info rows
anno_all <- anno_all[c(1:2187),] #Get rid of empty rows at the end
column_names <- c("exclusion","Type.Exp","WTA.plate","Look.Exp","original.384wp","LCM.exp","Pairs2","event2","LCM_found_conf2","key_pairs2","MN rupture timing (frame # to total # frames)",".1","capture time","division timing (division frame #, (total # frames))","Lysis timing (lysis time after time of last cell division)","notes1","Cells batch","Cell Cycle (\"h\" at lysis)","Pair","Chromosomal event (infered by analysis)","Notes on analysis","Seg CN Ratio","exclusion criteria","Genes detected >1","Genes detected >5","video","imaging.plate","WTA.plate2","SS2.Exp","lysis","protocol2","cDNA.amp2","notes2","BioA","cDNA.QC","Sequencing.A","Nextera_plate","Sequencing.B","SeqQC_total_reads","SeqQC_%_exonic_reads","SeqQC_#_genes","further.cDNA.amp","BioA (further cDNA amp)")
colnames(anno_all) <- column_names
setkey(anno_all, "WTA.plate")
anno_all <- anno_all[which(anno_all$WTA.plate != ""),]

#Now lets load the old full anno list to merge
anno_full <- data.table(read.csv(sprintf("%s/work_in_progress/Annotation_list_long_vnikos_1_9_21.csv", datadir)))
columns <- colnames(anno_full)[c(1:25,35:36,165:168)]
anno_full <- anno_full[,..columns] #Keep only the useful columns
setkey(anno_full, "WTA.plate")

#Setkeys for shared columns
cols <- intersect(colnames(anno_all), colnames(anno_full))
setkeyv(anno_all, cols)
setkeyv(anno_full, cols)

#merge files by WTA.plate name
anno_final <- merge(anno_all, anno_full, by = cols, all = T)

#Clean up the combined file
cols_ordered <- c("exclusion", "exclusion criteria", "WTA.plate", "WTA.plate2", "NTA.plate", "Nextera_plate", "Cells batch", "Pair", "Pairs", "Pairs2", "key_pairs", "key_pairs2", "control_group", "event", "event2", "comments_by_stam", "notes", "notes1", "notes2", "method.notes", "treatment", "Type.Exp", "LCM.exp", "Look.Exp", "LookRNAseq.Exp", "SS2.Exp", "experimental_label", "protocol", "protocol2", "LCM_found_conf", "LCM_found_conf2", "cDNA.amp", "cDNA.amp2", "cDNA.QC", "further.cDNA.amp", "BioA (further cDNA amp)", "BioA", "capture time", "lysis", "Lysis timing (lysis time after time of last cell division)", "Cell Cycle (\"h\" at lysis)", "MN rupture timing (frame # to total # frames)", "Chromosomal event (infered by analysis)", "Notes on analysis", "Seg CN Ratio", "imaging.plate", "video", "division timing (division frame #, (total # frames))", "Sequencing.A", "Sequencing.B", "SeqQC_total_reads", "SeqQC_%_exonic_reads", "SeqQC_#_genes", "mainGroup", "abbr_group", "cellType1", "cellType2", "cell.type", "cell.type.1", "original.384wp", "Orig.Imag.384wp", "final.384wp", "Seq.date", "Seq.platform", "Seq.kit", "filePath", "filePath_2nd", "Fastq_files", "Fastq_files2", "Fastq_files3") #reorder columns
anno_final <- anno_final[,..cols_ordered] #reorder the columns
na.fill(anno_final, fill = "")
anno_final[anno_final == ""] <- NA
anno_final$LCM_found_conf <- as.character(anno_final$LCM_found_conf)

#Add most recent QC  to the new anno list
all_QC <- data.table(readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", wkpath)))
setnames(all_QC, "id", "WTA.plate")
setkey(all_QC, "WTA.plate")
anno_final_QC <- merge(anno_final, all_QC, by = "WTA.plate", all = T)

write.csv(anno_final_QC, sprintf("%s/work_in_progress/Annotation_list_Oct13_2021.csv", datadir), row.names = FALSE)
anno_final_QC <- read.csv(sprintf("%s/work_in_progress/Annotation_list_Oct13_2021.csv", datadir))

#get the order of the original table back. 
anno_all <- data.table(read.csv("/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/SS4_list_Oct13_St.csv", header = F))

colnames(anno_all) <- as.character(anno_all[8,]) #use colnames row as headers
anno_all2 <- anno_all[,c(1:43)] #get rid of empty columns
anno_all2 <- anno_all2[-c(1:9),] #Get rid of empty info rows
anno_all2 <- anno_all2[c(1:2187),] #Get rid of empty rows at the end
column_names <- c("exclusion","Type.Exp","WTA.plate","Look.Exp","original.384wp","LCM.exp","Pairs2","event2","LCM_found_conf2","key_pairs2","MN rupture timing (frame # to total # frames)",".1","capture time","division timing (division frame #, (total # frames))","Lysis timing (lysis time after time of last cell division)","notes1","Cells batch","Cell Cycle (\"h\" at lysis)","Pair","Chromosomal event (infered by analysis)","Notes on analysis","Seg CN Ratio","exclusion criteria","Genes detected >1","Genes detected >5","video","imaging.plate","WTA.plate2","SS2.Exp","lysis","protocol2","cDNA.amp2","notes2","BioA","cDNA.QC","Sequencing.A","Nextera_plate","Sequencing.B","SeqQC_total_reads","SeqQC_%_exonic_reads","SeqQC_#_genes","further.cDNA.amp","BioA (further cDNA amp)")
colnames(anno_all2) <- column_names
anno_all2 <- anno_all2[which(anno_all2$WTA.plate != ""),]

#transpose the table and select columns in the order of the old file and transpose back. 
tmp <- t(anno_final_QC)
colnames(tmp) <- tmp[1,]
missing_samps <- setdiff(colnames(tmp), anno_all2$WTA.plate)
original_order_samps <- append(anno_all2$WTA.plate, missing_samps)
tmp2 <- tmp[,original_order_samps]
reordered_anno_final_QC <- data.table(t(tmp2))

#Save the reordered file
write.csv(reordered_anno_final_QC, sprintf("%s/work_in_progress/Annotation_list_Oct13_2021.csv", datadir), row.names = FALSE)

