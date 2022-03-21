
#Generate the directories based on experiments
dir_paths <- list()
for (i in experiments) {dir_paths <- cbind(dir_paths, sprintf("%s/%s/STAR", wkpath, i))}

#list all files in each directory folder
files <- c()
for (i in 1:length(dir_paths)){files <- append(files, paste(dir_paths[i], list.files(path = as.character(dir_paths[i]), pattern = ".ReadsPerGene.out.tab"), sep = "/"))}

#List sample names
files_samples <- files
for (i in 1:length(anno$Fastq_files)){
  res <- grep(pattern = sprintf("_%s_", anno$Fastq_files[i]), ignore.case = T, value=T, x = files)
  if (length(res) < 1) {res <- grep(pattern =  anno$Fastq_files[i], ignore.case = T, value=T, x = files)}
  files_samples[which(files == res[1])] <- anno$WTA.plate[i]
  files_samples[which(files == res[2])] <- sprintf("%s_L2", anno$WTA.plate[i]) #if there is a second file it is because there were two lanes
}
files_samples_v2 <- append(files_samples[which(files_samples %in% anno$WTA.plate)], files_samples[which(files_samples %in% sprintf("%s_L2", anno$WTA.plate))])
files_v2 <- append(files[which(files_samples %in% anno$WTA.plate)], files[which(files_samples %in% sprintf("%s_L2", anno$WTA.plate))])

qcmat <- do.call(cbind, lapply(files_v2,
                               function(f) {
                                 message("Reading: ", f)
                                 read.table(f, header = F, sep = "\t", stringsAsFactors = F, row.names = 1)[,1,drop=F] }))

colnames(qcmat) <- files_samples_v2

qc <- qcmat[grep("^N_", rownames(qcmat)), ]
counts <- qcmat[-grep("^N_", rownames(qcmat)), ]

qc <- rbind(qc, geneCounts=colSums(counts))
list(counts=counts, qc=qc)

#Divide by cols by col sums and multiply by 100 to get a percent
qc2 <- round(sweep(qc, 2, colSums(qc), FUN = "/")*100)
#Calculate th1, th5, th10. Genes covered by at least 1, 5, or 10 reads.
qcbygene <- rbind(colSums(counts>0), colSums(counts>4), colSums(counts>9))

qc3 <- t(rbind(qc2, qcbygene, totalReads=colSums(qc)))
qcs <- cbind(id = rownames(qc3), qc3)
colnames(qcs) <- c("id", rownames(qc), "th1", "th5", "th10", "totalReads")

#make analysis directory based on input file path
mk_analysis_dir <- sprintf("mkdir -p  %s/QC", wkpath)
system(mk_analysis_dir)

#Write QC information for this experiment to new file path
write.csv(x = qcs, file = sprintf("%s/QC/QC.csv", wkpath)) 





