#rawTPMtoCNpred

require(data.table)
require(matrixStats)

#read in data
geneRanges <- readRDS(sprintf("%s/geneRanges.rds", datadir))
rsemtpm <- readRDS(sprintf("%s/aggregated_results/all_experiments_rsemtpm.rds", dirpath))
controlSampleIDs <- readRDS( sprintf("%s/aggregated_results/controlSampleIDs.rds", dirpath))

#get intersecting genes
common_genes <- intersect(rownames(rsemtpm), geneRanges$id)
rsemtpm <- rsemtpm[common_genes,]

#reformat rownames as a column in rsemtpm
rsemtpm <- cbind(rownames(rsemtpm), rsemtpm)
colnames(rsemtpm)[1] <- "id"

#merge annotations and values by "id" column and reorder
rsemtpm <- merge(geneRanges, rsemtpm, by = "id")
rsemtpm <- data.table(rsemtpm[order(seqnames, start, end)])


tpm <- data.matrix(rsemtpm[,-c(1:6)])
#cap tpm values in the matrix
rsemtpm2 <- tpm
#maxExp <- 1000 #maximum tpm value from literature
#minExp <- 5 #minimum tpm value
#rsemtpm2[rsemtpm2 > maxExp] <- NA
#rsemtpm2[rsemtpm2 < minExp] <- NA

#for (i in 1:ncol(rsemtpm2)) {
#  outliers <- boxplot(rsemtpm2[,i], plot=FALSE)$out
#  index <- which(rsemtpm2[,i] %in% outliers)
#  rsemtpm2[index,i] <- NA
#}

#Check the number of values in the correct range
#message("non NAs values = ", sum(!is.na(rsemtpm2)))

#get rows/genes which have nonzero expression in any cell
wi_rows <- which(rowSums(rsemtpm2>5, na.rm = T) > 500)
#wi_cols <- names(which(colSums(rsemtpm2>5, na.rm = T)>2250))
#wi_cols <- append(colnames(rsemtpm)[c(1:6)], wi_cols)

#only keep genes with nonzero expression in >1 cell
#rsemtpm3 <- data.table( cbind(rsemtpm[,c(1:6)], rsemtpm2) )
#rsemtpm4 <- data.table(rsemtpm3[wi_rows,]) #..wi_cols

rsemtpm_red <- rsemtpm2[wi_rows,]
rsemtpm_anno <- data.table(rsemtpm[wi_rows,c(1:6)]) 


#use limma::voom to inverse variance weight RNA-seq counting data
library(BiocManager)
BiocManager::install("limma")
require(limma)
?limma
#norm_rsemtpm_red <- calcNormFactors(rsemtpm_red, method="RLE")
voom_tpm <- limma::voom(rsemtpm_red, normalize.method = "cyclicloess", plot=T)
#ivw_tpm <- (voom_tpm[["weights"]]*voom_tpm[["E"]]) / (rowSums(voom_tpm[["weights"]]))
ivw_tpm <- cbind(rsemtpm_anno, voom_tpm[["E"]])

hist(ivw_tpm[,c(450)])
hist(adt.default[,c("181013A_5H")])



#Return to just values for log transform
dm <- data.matrix(rsemtpm4[,-c(1:6)])
#dm_control <- data.matrix(rsemtpm[wi_rows,..controlSampleIDs])

#log the data.matrix
log_dm <- log2(dm)
#log_control <- log2(dm_control)

#subtract col means
log_dm2 <- scale(log_dm, scale = F)

#unlog values
log_dm3 <- 2^log_dm2

#subtract row means
means <- rowMeans(log_dm3, na.rm = T)
log_dm4 <- sweep(log_dm3, 1, means, "/")

#reattach annotations
rsemtpm5 <- data.table(cbind(rsemtpm4[,c(2)], log_dm4))

#Get avg value for chr
rsemtpm6 <- rsemtpm5 %>% 
              group_by(seqnames) %>%
                summarise_all(mean, na.rm = TRUE)

#Center at 2 instead 1
rsemtpm6[,-c(1)] <- 2*rsemtpm6[,-c(1)]

#Save the object
saveRDS(rsemtpm6, sprintf("%s/aggregated_results/normalized_rsemtpm_bychr.rds", dirpath))

boxplot(rsemtpm6[,c(2:10)], outline = F)
barchart(as.matrix(rsemtpm6[,c(2:3)]))

##### Pick up here, I was working on what to subtract from the rows.I used to subtract the row mean but that was skewed so i tried controls and then had to leave. 


