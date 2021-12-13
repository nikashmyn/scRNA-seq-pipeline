#Script used to finding genes with Differential expression

library(data.table)
library(readxl)
library(robustbase)
library(stats)
library(matrixStats)
library(ggplot2)
source('/pellmanlab/nikos/Stam_Etai_Scripts/scripts/fromRawTPMtoExprsRatio.R')
etaipath <- "/singlecellcenter/etai/ExperimentsData/RNA/ngs/CNV12.2/data"
nikospath <- "/pellmanlab/stam_niko/data/processed_bam/aggregated_results"
anno <- data.table(read_excel("/pellmanlab/nikos/Stam_Etai_Data/work_in_progress/Stamatis_list_v15_190910_edit.xlsx"))

controlSampleIDs2 <- readRDS(sprintf("%s/controlSampleIDs2.rds", nikospath))
controlSampleIDs <- readRDS(sprintf("%s/controlSampleIDs.rds", nikospath))

tpm_etai <- readRDS(sprintf("%s/rsemtpm.rds", etaipath))
tpm_nikos <- readRDS(sprintf("%s/all_experiments_rsemtpm.rds", nikospath))

Etai_Samples <- colnames(tpm_etai)
Nikos_Samples <- colnames(tpm_nikos)

Missing_Samples <- setdiff(Etai_Samples, Nikos_Samples)
All_Samples_v15 <- anno$WTA.plate

#Take note of missing samples
rows <- c()
for (i in 1:length(Missing_Samples)) {
  for (j in 1:length(All_Samples_v15)) {
    if (Missing_Samples[i]==All_Samples_v15[j]) {
      rows <- append(rows, j)
    }
  }
}
anno_missing_samples <- anno[rows, ]

#Intersect and reduce etai's samples with nikos'
columns <- c()
for (i in 1:length(Missing_Samples)) {
  for (j in 1:length(colnames(tpm_etai))) {
    if (Missing_Samples[i]==colnames(tpm_etai)[j]) {
      columns <- append(columns, j)
    }
  }
}
tpm_etai_new <- tpm_etai[,-columns]
tpm_etai_new <- tpm_etai_new[,order(colnames(tpm_etai_new))]

#Intersect and reduce nikos samples with etai's
columns <- c()
for (i in 1:length(colnames(tpm_etai_new))) {
  for (j in 1:length(colnames(tpm_nikos))) {
    if (colnames(tpm_etai_new)[i]==colnames(tpm_nikos)[j]) {
      columns <- append(columns, j)
    }
  }
}
tpm_nikos_new <- tpm_nikos[,columns]
tpm_nikos_new <- tpm_nikos_new[,order(colnames(tpm_nikos_new))]

#Intersect and reduce etai's genes with nikos'
rows_missing <- setdiff( rownames(tpm_etai_new), rownames(tpm_nikos_new) )
rows <- c()
for (i in 1:length(rows_missing)) {
  for (j in 1:length(rownames(tpm_etai_new))) {
    if (rows_missing[i]==rownames(tpm_etai_new)[j]) {
      rows <- append(rows, j)
    }
  }
}
tpm_etai_new <- tpm_etai_new[-rows,]

#Exclude lowly expressed genes based on control samples
minDetectionLevel = 5 #option1: 5 #original: 5
minNumOfSamplesToDetect = 3 #option1: 10 #original: 3
wi <- which(rowSums(tpm_nikos_new[, controlSampleIDs] >= minDetectionLevel) >= minNumOfSamplesToDetect)
tpm_nikos_new <- tpm_nikos_new[wi,]

#import geneRanges annotations
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", nikospath))
genes <- geneRanges$id

#Intersect and reduce geneRanges and tpm genes
final_genes <- intersect(rownames(tpm_nikos_new), genes)
tpm_nikos_new <- tpm_nikos_new[final_genes,]
tpm_etai_new <- tpm_etai_new[final_genes,]
geneRanges <- geneRanges[geneRanges$id %in% final_genes]

#Final tpm dataframe
rsemtpm_nikos <- cbind(geneRanges, tpm_nikos_new)
rsemtpm_etai <- cbind(geneRanges, tpm_etai_new)

rsemtpm_nikos <- readRDS("/pellmanlab/stam_niko/data/comparable_datasets/nikos_dataset/aggregated_results2/adt_nolog_nocent.rds")

#rsemtpm_nikos <- cbind(rsemtpm_nikos[,c(1:4)], rsemtpm_nikos[,..controlSampleIDs])

#### graph 1 #####
#Choose the gene
gene <- ""
gene_index <- 2

#Get the sample distribution
rsemtpm_nikos_chr1 <- rsemtpm_nikos[rsemtpm_nikos$seqnames == 'chr1']
rsemtpm_nikos_geneordered <- rsemtpm_nikos_chr1[order(rsemtpm_nikos_chr1$start),]
samples_dist <- rsemtpm_nikos_geneordered[gene_index,-c(1:4)]
samples_dist <- as.numeric(samples_dist)

#Get the max tpm (plus one so I can bin up to this amount and include the max sample)
max_tpm <- max(samples_dist) + 1

#Create X bins from the max tpm
number_of_bins <- 100
size_of_bin <- max_tpm / number_of_bins

#Create the bin divisions
bin_divs <- c(0)
counter <- 0
for (i in 1:number_of_bins) {
  counter <- counter + size_of_bin
  bin_divs <- append(bin_divs, counter)
}

#Count the number of samples in each tpm bin
samples_dist_binned_freq <- c()
for ( i in 2:length(bin_divs) ) {
  tmp <- samples_dist[samples_dist > bin_divs[i-1]]
  tmp2 <- tmp[tmp < bin_divs[i]]
  freq <- length(tmp2)
  samples_dist_binned_freq <- append(samples_dist_binned_freq, freq)
}

#plot
barplot(samples_dist_binned_freq,
        main="Sample Distribution of TPM for a Gene", 
        xlab=sprintf("TPM (%s bins of size %s)", number_of_bins, size_of_bin),
        ylab="Frequency (Samples)")


### Graph 2 ### 
number_of_bins <- 20

#Get the sample distribution
rsemtpm_nikos_chr1 <- rsemtpm_nikos[rsemtpm_nikos$seqnames == 'chr1']
rsemtpm_nikos_geneordered <- rsemtpm_nikos_chr1[order(rsemtpm_nikos_chr1$start),]
rsemtpm_nikos_ids <- rsemtpm_nikos_geneordered$id

freq_table <- c(1:number_of_bins)
for ( i in 1:nrow(rsemtpm_nikos_geneordered)) {
  
  rowname <- rownames(rsemtpm_nikos_geneordered)
  samples_dist <- rsemtpm_nikos_geneordered[i,-c(1:4)]
  samples_dist <- as.numeric(samples_dist)
  
  #Get the max tpm (plus one so I can bin up to this amount and include the max sample)
  max_tpm <- max(samples_dist) + 1
  
  #Create X bins from the max tpm
  size_of_bin <- max_tpm / number_of_bins
  
  #Create the bin divisions
  bin_divs <- c(1)
  counter <- 0
  for (i in 1:number_of_bins) {
    counter <- counter + size_of_bin
    bin_divs <- append(bin_divs, counter)
  }
  
  #Count the number of samples in each tpm bin
  samples_dist_binned_cnts <- c()
  for ( i in 2:length(bin_divs) ) {
    tmp <- samples_dist[samples_dist > bin_divs[i-1]]
    tmp2 <- tmp[tmp < bin_divs[i]]
    counts <- length(tmp2)
    samples_dist_binned_cnts <- append(samples_dist_binned_cnts, counts)
  }
  
  #Normalize counts to frequency
  xmax <- max(samples_dist_binned_cnts)
  xmin <- min(samples_dist_binned_cnts)
  samples_dist_binned_freq <- (samples_dist_binned_cnts - xmin) / (xmax - xmin)
  freq_table <- rbind(freq_table, samples_dist_binned_freq)
}

colnames(freq_table) <- bin_divs[-c(1)]
freq_table <- freq_table[-c(1),]
rownames(freq_table) <- rsemtpm_nikos_ids

heatmap(freq_table[,], Colv = NA, #Rowv = NA,
        labRow = NA, #labCol = NA,
        main = "TPM Frequency Heatmap", 
        xlab = "TPM", 
        ylab = "Genes in Chr1")

row.clusters = hclust(dist(freq_table)) # get clusters
clusters <- cutree(row.clusters,k=3) # break into k=3 clusters
freq_table_clustered <- data.table(cbind(freq_table, clusters), keep.rownames = "id") 
freq_table_cluster <- freq_table_clustered[(clusters == 2 | clusters == 3)] #get third cluster
gene_cluster_ids <- freq_table_cluster[,c(1)]  #set rownames as genes
freq_table_cluster_matrix <- as.matrix(freq_table_cluster[,-c("clusters")], rownames = 1) #remove id col

saveRDS(gene_cluster_ids, "/pellmanlab/stam_niko/data/comparable_datasets/nikos_dataset/aggregated_results2/differentially_expr_genes.rds")


heatmap(freq_table_cluster_matrix[,], Colv = NA, #Rowv = NA,
        labRow = NA, #labCol = NA, 
        main = "TPM Frequency Heatmap", 
        xlab = "TPM", 
        ylab = "Genes in Chr1")

#freq_kmeans <- kmeans(freq_table, 2)


### graph 3 ### Get clustering by DE in this part
number_of_bins <- 100

#Get the sample distribution
rsemtpm_nikos_chr1 <- rsemtpm_nikos[rsemtpm_nikos$seqnames == 'chr1',]
rsemtpm_nikos_geneordered <- rsemtpm_nikos_chr1[order(rsemtpm_nikos_chr1$start),]
rsemtpm_nikos_geneordered <- rsemtpm_nikos_geneordered[id %in% gene_cluster_ids$id,]

rsemtpm_nikos_geneordered <- data.matrix(data.frame(sapply(rsemtpm_nikos_geneordered[,-c(1:4)], function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}})))

freq_table <- c(1:number_of_bins)
for ( i in 1:ncol(rsemtpm_nikos_geneordered) ) {
  genes_dist <- rsemtpm_nikos_geneordered[,i]
  genes_dist <- as.numeric(genes_dist)
  
  #Get the max tpm (plus one so I can bin up to this amount and include the max sample)
  max_tpm <- max(genes_dist) + 1
  
  #Create X bins from the max tpm
  size_of_bin <- max_tpm / number_of_bins
  
  #Create the bin divisions
  bin_divs <- c(1)
  counter <- 0
  for (i in 1:number_of_bins) {
    counter <- counter + size_of_bin
    bin_divs <- append(bin_divs, counter)
  }
  
  #Count the number of samples in each tpm bin
  genes_dist_binned_cnts <- c()
  for ( i in 2:length(bin_divs) ) {
    tmp <- genes_dist[genes_dist > bin_divs[i-1]]
    tmp2 <- tmp[tmp < bin_divs[i]]
    counts <- length(tmp2)
    genes_dist_binned_cnts <- append(genes_dist_binned_cnts, counts)
  }
  
  #Normalize counts to frequency
  xmax <- max(genes_dist_binned_cnts)
  xmin <- min(genes_dist_binned_cnts)
  genes_dist_binned_freq <- (genes_dist_binned_cnts - xmin) / (xmax - xmin)
  freq_table <- rbind(freq_table, genes_dist_binned_freq)
}

colnames(freq_table) <- freq_table[1,]
freq_table <- freq_table[-c(1),]
rownames(freq_table) <- colnames(rsemtpm_nikos_geneordered)
freq_table <- freq_table[-c(212,474,480,482),]

heatmap(freq_table[,c(1:(number_of_bins/5))], Colv = NA,# Rowv = NA,
        #labRow = NA, labCol = NA,
        main = "Gene Frequency Heatmap", 
        xlab = "TPM", 
        ylab = "Cells")

#freq_table <- freq_table[-c("X180814A_5D"),]

row.clusters = hclust(dist(freq_table)) # get clusters
clusters <- cutree(row.clusters,k=6) # break into k=3 clusters
freq_table_clustered <- data.table(cbind(freq_table, clusters), keep.rownames = "id") 
freq_table_cluster <- freq_table_clustered[clusters == 1] # 1-3 are good clusters
freq_table_cluster_ids <- freq_table_cluster[,c(1)]  #set rownames as genes
freq_table_cluster_matrix <- as.matrix(freq_table_cluster[,-c("clusters")], rownames = 1) #remove id col

heatmap(freq_table_cluster_matrix[,c(1:(number_of_bins/5))], Colv = NA, #Rowv = NA,
        #labRow = NA, labCol = NA,
        main = "TPM Frequency Heatmap", 
        xlab = "TPM", 
        ylab = "Genes in Chr1")