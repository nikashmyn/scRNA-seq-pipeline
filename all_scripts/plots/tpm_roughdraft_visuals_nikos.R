#rough draft visual ideas for raw tpm data

library(data.table)
library(readxl)
library(robustbase)
library(stats)
library(matrixStats)
source('/pellmanlab/nikos/Stam_Etai_Scripts/fromRawTPMtoExprsRatio.R')
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

adt_etai <- readRDS(sprintf("%s/adt.rds", etaipath))
adt_nikos <- readRDS(sprintf("%s/adt.rds", nikospath))
#adt.zscore.rds seems wrong it is all too big and not distributed around zero. 

adt_etai <- as.data.frame(adt_etai)
adt_nikos <- as.data.frame(adt_nikos)

rownames(adt_etai) <- adt_etai[,1]
rownames(adt_nikos) <- adt_nikos[,1]

adt_nikos_short <- adt_nikos[rownames(adt_etai) ,]
adt_etai_short <- adt_etai

shared_columns <- intersect(colnames(adt_etai_short), colnames(adt_nikos_short))

adt_nikos_short <- adt_nikos_short[, shared_columns]
adt_etai_short <- adt_etai_short[, shared_columns]

adt_nikos_short_noanno <- adt_nikos_short[,-c(1:4)]
adt_etai_short_noanno <- adt_etai_short[,-c(1:4)]

adt_nikos_short_ordered <- adt_nikos_short_noanno[,order(colnames(adt_nikos_short_noanno))]
adt_etai_short_ordered <- adt_etai_short_noanno[,order(colnames(adt_etai_short_noanno))]

adt_nikos_short_ordered_wanno <- cbind(adt_nikos_short[,c(1:4)], adt_nikos_short_ordered)
adt_etai_short_ordered_wanno <- cbind(adt_etai_short[,c(1:4)], adt_etai_short_ordered)

adt_nikos_final <- adt_nikos_short_ordered_wanno[order(rownames(adt_nikos_short_ordered_wanno)),]
adt_etai_final <- adt_etai_short_ordered_wanno[order(rownames(adt_etai_short_ordered_wanno)),]

#------------------------------------------------------------------------------
#Graph 0.1: Raw Signal
#This is the raw gene expression for chr1 with gene quantity based bins.

#NOTES: 
#1) For chr1 etai has 19 bins and he says there are 50 genes each. That cannot be the case. There are #5194 genes in chr1. 5194genes/19bins = 273.4 genes/bin There must be more complexity like genes with #zero expression in all samples are thrown out or something.
#2) I believe etai gets rid of some genes. I need to figure out how he does this. What are the criteria?

rsemtpm_nikos_chr1 <- rsemtpm_nikos[rsemtpm_nikos$seqnames == 'chr1']
#look at this
#nonzeros <- nonzeros[order(seqnames, start)]

rsemtpm_nikos_geneordered <- rsemtpm_nikos_chr1[order(rsemtpm_nikos_chr1$start),]


number_of_genes <- nrow(rsemtpm_nikos_geneordered)
number_of_bins <- 20

bin_size <- number_of_genes/number_of_bins

genome_bins <- as.integer(ceiling((number_of_genes) / bin_size))

bins <- c(0)
counter <- 0
for (i in 1:genome_bins) {
  counter <- counter + bin_size
  bins <- append(bins, counter)
}

rsemtpm_nikos_chr1_binsummed <- colnames(rsemtpm_nikos_geneordered[,-c(1:6)])
for ( i in 2:length(bins) ) {
  tmp <- rsemtpm_nikos_geneordered[bins[i-1]:bins[i],]
  col_sums <- colSums(tmp[,-c(1:6)])
  rsemtpm_nikos_chr1_binsummed <- rbind(rsemtpm_nikos_chr1_binsummed, col_sums)
}
rsemtpm_nikos_chr1_binsummed <- rsemtpm_nikos_chr1_binsummed[-c(1),]
rownames(rsemtpm_nikos_chr1_binsummed) <- NULL
rsemtpm_nikos_chr1_binsummed <- as.data.frame(rsemtpm_nikos_chr1_binsummed)

df <- data.frame(sapply(rsemtpm_nikos_chr1_binsummed, function(x) { if(is.numeric(x) == F) {as.numeric(as.character(x))} else {x}}))

df2 <- t(df)

boxplot(df2, use.cols = TRUE, outline=FALSE,
        main="Raw Gene Expression for Chr1", 
        xlab=sprintf("Gene Bins (%s genes/bin)", as.integer(bin_size)), ylab="Expression (TPM)")

points(df2[15,], col = 150)

#--------------------------------------------------------------------------------------
#Graph 1: z-score boxplots with with sample points

#161228_A2 monosomy chr1 a control sample
#170208_A1 disomy control sample
#181012_5B disomy control sample
#190628_5A
outline_switch = F
sample <- "161228_A2"
chromosome <- 'chr1'
number_of_bins <- 20

for (i in 1:2) {
  #get chr and order TPM subobject 
  rsemtpm_nikos_chr1 <- rsemtpm_nikos[rsemtpm_nikos$seqnames == chromosome]
  rsemtpm_nikos_geneordered <- rsemtpm_nikos_chr1[order(rsemtpm_nikos_chr1$start),]
  sample_index <- grep(sprintf("^%s$", sample), colnames(rsemtpm_nikos_geneordered[,-c(1:6)]))
  
  #Get the control sample indices
  rsemtpm_nikos_control_indices <- c()
  for (i in 1:ncol(rsemtpm_nikos_geneordered)) {
    if (colnames(rsemtpm_nikos_geneordered)[i] %in% shared_columns ) {
      rsemtpm_nikos_control_indices <- append(rsemtpm_nikos_control_indices, i)
    }
  }
  
  #if (control_switch == T) {
  #  #control sample TPMs object
  #  shared_columns <- intersect(colnames(rsemtpm_nikos_geneordered), controlSampleIDs)
  #  rsemtpm_nikos_geneordered <- as.data.frame(rsemtpm_nikos_geneordered, col.names = T)
  #  rsemtpm_nikos_control <- rsemtpm_nikos_geneordered[,shared_columns]
  #
  #  #Make all values into numerical values. Data manipulation affected consistency. 
  #  rsemtpm_nikos_control <- #data.matrix(data.frame(sapply(rsemtpm_nikos_control[,-c(1:6)],function(x) { #if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}})))
  #  
  #  #rsemtpm_nikos_control <- log2(rsemtpm_nikos_control + 1)
  #}
  
  rsemtpm_nikos_geneordered <- data.matrix(data.frame(sapply(rsemtpm_nikos_geneordered[,-c(1:6)], function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}})))
  
  #rsemtpm_nikos_geneordered <- log2(rsemtpm_nikos_geneordered + 1)
  
  #get row statistics
  if (control_switch == T) {
    means <- rowMeans(rsemtpm_nikos_control)
    sds <- rowSds(rsemtpm_nikos_control)
  } else {
    means <- rowMeans(rsemtpm_nikos_geneordered)
    sds <- rowSds(rsemtpm_nikos_geneordered)
  }
  
  #Replace all entries in df with zscore of that entry
  for ( i in 1:nrow(rsemtpm_nikos_geneordered )) {
    for ( j in 1:ncol(rsemtpm_nikos_geneordered) ) {
      rsemtpm_nikos_geneordered[i,j] <- (rsemtpm_nikos_geneordered[i,j] - means[i]) / sds[i]
    }
  }
  
  number_of_genes <- nrow(rsemtpm_nikos_geneordered)
  
  bin_size <- number_of_genes/number_of_bins
  
  genome_bins <- as.integer(ceiling((number_of_genes) / bin_size))
  
  bins <- c(0)
  counter <- 0
  for (i in 1:genome_bins) {
    counter <- counter + bin_size
    bins <- append(bins, counter)
  }
  
  rsemtpm_nikos_chr1_binmedian <- colnames(rsemtpm_nikos_geneordered)
  for ( i in 2:length(bins) ) {
    tmp <- rsemtpm_nikos_geneordered[bins[i-1]:bins[i],]
    abstmp <- abs(tmp)
    row_abssums <- as.vector(rowSums(abstmp, na.rm = T))
    whichmedian <- function(x) which.min(abs(x - median(x)))
    median_row_index <- whichmedian(row_abssums)
    #median_row_index <- which.min(row_abssums)
    #median_row_index <- which.max(row_abssums)
    #col_med <- tmp[median_row_index,]
    #tmp_medial_genes <- tmp[2-median_row_index:median_row_index+2,]
    #col_med <- colMedians(tmp_medial_genes)
    col_med <- colMedians(tmp)
    #col_med <- colMeans(tmp)
    rsemtpm_nikos_chr1_binmedian <- rbind(rsemtpm_nikos_chr1_binmedian, col_med)
  }
  rsemtpm_nikos_chr1_binmedian <- rsemtpm_nikos_chr1_binmedian[-c(1),]
  rownames(rsemtpm_nikos_chr1_binmedian) <- NULL
  rsemtpm_nikos_chr1_binmedian <- as.data.frame(rsemtpm_nikos_chr1_binmedian)
  
  df <- data.frame(sapply(rsemtpm_nikos_chr1_binmedian, function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}}))
  
  df2 <- t(df)
  
  if (outline_switch == T) {
    title = "Standardized Expression for Chr1 (w/ Outliers)"
    ytitle = "Z-score of Medially Variant Gene in Bin (Stds)"
  } else {
    title = "Standardized Expression for Chr1 "
    ytitle = "Z-score of Medially Variant Gene in Bin (Stds)"
  }
  
  # replace df2 with df2[rsemtpm_nikos_control_indices,] control boxplots
  boxplot(df2[rsemtpm_nikos_control_indices,], use.cols = TRUE, outline=outline_switch,
          main=title, 
          xlab=sprintf("Gene Bins (%s genes/bin)", as.integer(bin_size)), 
          ylab=ytitle)
  #ylim=c(-.7,.7))
  grid(nx = NULL, ny = NULL, col = "lightgray",
       lwd = par("lwd"), equilogs = TRUE)
  points(df2[sample_index,], col = 150)
  
  #Switch outline on for second loop
  outline_switch = T
}

#------------------------------------------------------------------------------
#Graph 2: Fractional Scatter Plot with quartiles
##NOTES:
#How do I deal with the zero row medians when I get the fraction. Means work because they are usually #non-zero. But medians have many zeros. 

#161228_A2 monosomy chr1 a control sample
#170208_A1 disomy control sample
#181012_5B
#190628_5A
control_switch = F
sample <- "161228_A2"
chromosome <- 'chr1'
number_of_bins <- 20

for (i in 1:2) {
  #get chr and order TPM subobject 
  rsemtpm_nikos_chr1 <- rsemtpm_nikos[rsemtpm_nikos$seqnames == chromosome] # turn back to rsemtpm_nikos
  rsemtpm_nikos_geneordered <- rsemtpm_nikos_chr1[order(rsemtpm_nikos_chr1$start),]
  sample_index <- grep(sprintf("^%s$", sample), colnames(rsemtpm_nikos_geneordered[,-c(1:6)]))
  
  if (control_switch == T) {
    #control sample TPMs object
    shared_columns <- intersect(colnames(rsemtpm_nikos_geneordered), controlSampleIDs)
    rsemtpm_nikos_geneordered <- as.data.frame(rsemtpm_nikos_geneordered, col.names = T)
    rsemtpm_nikos_control <- rsemtpm_nikos_geneordered[,shared_columns]
    
    #Make all values into numerical values. Data manipulation affected consistency. 
    rsemtpm_nikos_control <- data.matrix(data.frame(sapply(rsemtpm_nikos_control[,-c(1:6)],function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}})))
    
    #Lets make all the values log scale to reduce the effect of the outliers
    #rsemtpm_nikos_control <- log2(rsemtpm_nikos_control + 1)
  }
  
  #Make all values into numerical values. Data manipulation affected consistency. 
  rsemtpm_nikos_geneordered <- data.matrix(data.frame(sapply(rsemtpm_nikos_geneordered[,-c(1:6)], function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}})))
  
  #Lets make all the values log scale to reduce the effect of the outliers
  #rsemtpm_nikos_geneordered <- log2(rsemtpm_nikos_geneordered + 1)
  
  if (control_switch == F) {meds <- rowMedians(rsemtpm_nikos_geneordered)} else {meds <- rowMedians(rsemtpm_nikos_control)}
  
  for ( i in 1:nrow(rsemtpm_nikos_geneordered )) {
    if (meds[i] == 0) {
      meds[i] <- NaN
    }
    for ( j in 1:ncol(rsemtpm_nikos_geneordered) ) {
      rsemtpm_nikos_geneordered[i,j] <- rsemtpm_nikos_geneordered[i,j]/meds[i]
    }
  }
  
  number_of_genes <- nrow(rsemtpm_nikos_geneordered)
  
  bin_size <- number_of_genes/number_of_bins
  
  genome_bins <- as.integer(ceiling((number_of_genes) / bin_size))
  
  bins <- c(0)
  counter <- 0
  for (i in 1:genome_bins) {
    counter <- counter + bin_size
    bins <- append(bins, counter)
  }
  
  rsemtpm_nikos_chr1_binmedian <- colnames(rsemtpm_nikos_chr1[,-c(1:6)])
  rsemtpm_nikos_chr1_binq2 <- colnames(rsemtpm_nikos_chr1[,-c(1:6)])
  rsemtpm_nikos_chr1_binq4 <- colnames(rsemtpm_nikos_chr1[,-c(1:6)])
  for ( i in 2:length(bins) ) {
    tmp <- rsemtpm_nikos_geneordered[bins[i-1]:bins[i],]
    q2 <- colQuantiles(tmp, na.rm = T)[,c(2)]
    #col_med <- colMedians(tmp, na.rm = T)
    col_med <- colQuantiles(tmp, na.rm = T)[,c(3)]
    q4 <- colQuantiles(tmp, na.rm = T)[,c(4)]
    rsemtpm_nikos_chr1_binmedian <- rbind(rsemtpm_nikos_chr1_binmedian, col_med)
    rsemtpm_nikos_chr1_binq2 <- rbind(rsemtpm_nikos_chr1_binq2, q2)
    rsemtpm_nikos_chr1_binq4 <- rbind(rsemtpm_nikos_chr1_binq4, q4)
  }
  rsemtpm_nikos_chr1_binmedian <- rsemtpm_nikos_chr1_binmedian[-c(1),]
  rsemtpm_nikos_chr1_binq2 <- rsemtpm_nikos_chr1_binq2[-c(1),]
  rsemtpm_nikos_chr1_binq4 <- rsemtpm_nikos_chr1_binq4[-c(1),]
  rownames(rsemtpm_nikos_chr1_binmedian) <- NULL
  rownames(rsemtpm_nikos_chr1_binq2) <- NULL
  rownames(rsemtpm_nikos_chr1_binq4) <- NULL 
  rsemtpm_nikos_chr1_binmedian <- as.data.frame(rsemtpm_nikos_chr1_binmedian)
  rsemtpm_nikos_chr1_binq2 <- as.data.frame(rsemtpm_nikos_chr1_binq2)
  rsemtpm_nikos_chr1_binq4 <- as.data.frame(rsemtpm_nikos_chr1_binq4)
  df <- data.frame(sapply(rsemtpm_nikos_chr1_binmedian, function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}}))
  df.q2 <- data.frame(sapply(rsemtpm_nikos_chr1_binq2, function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}}))
  df.q4 <- data.frame(sapply(rsemtpm_nikos_chr1_binq4, function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}}))
  
  df2 <- t(df)
  df2.q2 <- t(df.q2)
  df2.q4 <- t(df.q4)
  
  if (control_switch == T) {
    title = "(Control) Fractional Expression Distribution for Chr1"
    ytitle = "Raw TPM / Control Median TPM (TPM/TPM)"
  } else {
    title = "Fractional Expression Distribution for Chr1"
    ytitle = "Raw TPM / Median TPM (TPM/TPM)"
  }
  
  plot(df2[sample_index,], col = 150,
       ylim=c(-.1, 2.5),
       main=title, 
       xlab=sprintf("Gene Bins (%s genes/bin)", as.integer(bin_size)),
       ylab=ytitle)
  abline(a=1, b=0, lty="dashed")
  grid(nx = NULL, ny = NULL, col = "lightgray",
       lwd = par("lwd"), equilogs = TRUE)
  arrows(x0=c(1:20), y0=df2.q2[sample_index,], x1=c(1:20), y1=df2.q4[sample_index,], code=3, angle=90, length=0.05, col="blue", lwd=2)
  
  #Switch control based median on for second loop
  control_switch = T
}

#Graph 3: Normalization

sample <- "170208_A1"
chromosome <- 'chr1'
number_of_bins <- 20

#get chr and order TPM subobject 
rsemtpm_nikos_chr1 <- rsemtpm_nikos[rsemtpm_nikos$seqnames == chromosome] # turn back to rsemtpm_nikos
rsemtpm_nikos_geneordered <- rsemtpm_nikos_chr1[order(rsemtpm_nikos_chr1$start),]
sample_index <- grep(sprintf("^%s$", sample), colnames(rsemtpm_nikos_geneordered[,-c(1:6)]))

#Make all values into numerical values. Data manipulation affected consistency. 
rsemtpm_nikos_geneordered <- data.matrix(data.frame(sapply(rsemtpm_nikos_geneordered[,-c(1:6)], function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}})))

#Lets make all the values log scale to reduce the effect of the outliers
rsemtpm_nikos_geneordered <- log2(rsemtpm_nikos_geneordered + 1)

#Calculate Statistics for normalization
max <- rowMaxs(rsemtpm_nikos_geneordered)
min <- rowMins(rsemtpm_nikos_geneordered)

#Perform max-min scaling
for ( i in 1:nrow(rsemtpm_nikos_geneordered )) {
  for ( j in 1:ncol(rsemtpm_nikos_geneordered) ) {
    rsemtpm_nikos_geneordered[i,j] <- (rsemtpm_nikos_geneordered[i,j] - min[i]) / (max[i] - min[i])
  }
}

#Create bins based on dataframe size
number_of_genes <- nrow(rsemtpm_nikos_geneordered)
bin_size <- number_of_genes/number_of_bins
genome_bins <- as.integer(ceiling((number_of_genes) / bin_size))
bins <- c(0)
counter <- 0
for (i in 1:genome_bins) {
  counter <- counter + bin_size
  bins <- append(bins, counter)
}

#Bin the dataframe and get bin statistics
rsemtpm_nikos_chr1_binmedian <- colnames(rsemtpm_nikos_chr1[,-c(1:6)])
rsemtpm_nikos_chr1_binq2 <- colnames(rsemtpm_nikos_chr1[,-c(1:6)])
rsemtpm_nikos_chr1_binq4 <- colnames(rsemtpm_nikos_chr1[,-c(1:6)])
for ( i in 2:length(bins) ) {
  tmp <- rsemtpm_nikos_geneordered[bins[i-1]:bins[i],]
  q2 <- colQuantiles(tmp, na.rm = T)[,c(2)]
  col_med <- colQuantiles(tmp, na.rm = T)[,c(3)]
  q4 <- colQuantiles(tmp, na.rm = T)[,c(4)]
  rsemtpm_nikos_chr1_binmedian <- rbind(rsemtpm_nikos_chr1_binmedian, col_med)
  rsemtpm_nikos_chr1_binq2 <- rbind(rsemtpm_nikos_chr1_binq2, q2)
  rsemtpm_nikos_chr1_binq4 <- rbind(rsemtpm_nikos_chr1_binq4, q4)
}

#Reformat the new binned dataframe
rsemtpm_nikos_chr1_binmedian <- rsemtpm_nikos_chr1_binmedian[-c(1),]
rsemtpm_nikos_chr1_binq2 <- rsemtpm_nikos_chr1_binq2[-c(1),]
rsemtpm_nikos_chr1_binq4 <- rsemtpm_nikos_chr1_binq4[-c(1),]
rownames(rsemtpm_nikos_chr1_binmedian) <- NULL
rownames(rsemtpm_nikos_chr1_binq2) <- NULL
rownames(rsemtpm_nikos_chr1_binq4) <- NULL 
rsemtpm_nikos_chr1_binmedian <- as.data.frame(rsemtpm_nikos_chr1_binmedian)
rsemtpm_nikos_chr1_binq2 <- as.data.frame(rsemtpm_nikos_chr1_binq2)
rsemtpm_nikos_chr1_binq4 <- as.data.frame(rsemtpm_nikos_chr1_binq4)
df <- data.frame(sapply(rsemtpm_nikos_chr1_binmedian, function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}}))
df.q2 <- data.frame(sapply(rsemtpm_nikos_chr1_binq2, function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}}))
df.q4 <- data.frame(sapply(rsemtpm_nikos_chr1_binq4, function(x) { if(is.numeric(x)==F) {as.numeric(as.character(x))} else {x}}))


df2 <- t(df)
df2.q2 <- t(df.q2)
df2.q4 <- t(df.q4)

title = "Max-Min Normalized log2(TPM) Distribution for Chr1"
ytitle = "Normlized (0-1) log2(TPM) [Unitles]"

plot(df2[sample_index,], col = 150,
     ylim=c(-.001, 1),
     main=title, 
     xlab=sprintf("Gene Bins (%s genes/bin)", as.integer(bin_size)),
     ylab=ytitle)
arrows(x0=c(1:20), y0=df2.q2[sample_index,], x1=c(1:20), y1=df2.q4[sample_index,], code=3, angle=90, length=0.05, col="blue", lwd=2)