#########################################################
### Contour Map Function (TPM Ratio vs Variant Ratio) ###
#########################################################

#vector difference between alleles cnts per chr
saveRDS(allele_diff_mat.bychr, file=sprintf("%s/aggregated_results/allele_diff_mat.bychr.rds", dirpath))
#absolute difference between allele cnts per chr
saveRDS(abs_allele_diff_mat.bychr, file=sprintf("%s/aggregated_results/abs_allele_diff_mat.bychr.rds", dirpath))
#total SNP counts over both alleles normalized together after summing cnts. 
saveRDS(var_mat.bychr, file=sprintf("%s/aggregated_results/var_mat.bychr.rds", dirpath))
#total SNP counts over both alleles taken after separate allele normalization. ~same as var_mat.bychr. 
saveRDS(var_mat_sep.bychr, file=sprintf("%s/aggregated_results/var_mat_sep.bychr.rds", dirpath))

#QC Data Frame and high QC IDs based on 5 read cnt threshold for genes
all_QC <- readRDS(file = sprintf("%s/aggregated_results/all_QC.rds", dirpath))
high_qc_ids <- as.character(all_QC[which(all_QC$th5 >= quantile(all_QC$th5, c(.10))),]$id) #Excludes bottom 10%

###################
### Contour Map ###
###################

#get families from anno list
setkey(anno, WTA.plate)
families <- sort(table(anno[Pairs != "NA" & Pairs != "mother"]$Pairs), decreasing = T)

#Set seed for reproducability
set.seed(123)

#Prepare Visual Objects
myfamily = names(families)[142] #43, 24, 142, 30, 25
message(myfamily)
myids <- anno[Pairs %in% myfamily]$WTA.plate
visual.data <- data.table(TPM = as.numeric(flatten(adt.bychr[,..myids])))
visual.data$VAR <- as.numeric(flatten(var_mat.bychr[,..myids]))
visual.data$DIFF <- as.numeric(flatten(varabs_mat.bychr[,..myids]))
visual.data$chr <- as.factor(rep(unique(adt.bychr$seqnames), length = length(visual.data$VAR)))
visual.data$cell <- rep(myids, each=23-length(exclude_chrs))

####################
### Cluster Plot ###
####################
#shorten data name
df <- visual.data[, c("TPM", "VAR", "DIFF")]

#Visualize optimal number of clusters
wss <- fviz_nbclust(df, kmeans, method = "wss")
silhouette <- fviz_nbclust(df, kmeans, method = "silhouette")

#Extract optimal number of means
optimal_means_data <- fviz_nbclust(df, kmeans, method = "silhouette")$data
optimal_means <-as.numeric(optimal_means_data$clusters[which.max(optimal_means_data$y)])

#Visualize cluster map using optimal number of means
#cluster.data <- kmeans(df, centers = optimal_means, nstart = 25)
cluster.data <- kmeans(df, centers = 3, nstart = 25)

####################
### Final Visual ###
####################

#TPM/VAR plane plot
VAR.visual <- ggplot(data = visual.data, mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
  geom_point(mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
  xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "TPM Ratio", y = "SNP Count Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily))

#TPM/Diff plane plot
DIFF.visual <- ggplot(data = visual.data, mapping = aes(x = TPM, y = DIFF, color = chr, shape = cell)) +
  geom_point(mapping = aes(x = TPM, y = DIFF, color = chr, shape = cell)) +
  #xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "TPM Ratio", y = "Allelic Fraction Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily))

#Cluster visual
cluster.visual <- fviz_cluster(cluster.data, data = df, geom = "point")

grid.arrange(arrangeGrob(VAR.visual), arrangeGrob(cluster.visual, DIFF.visual, ncol=1, nrow=2), widths=c(2,1))







#########################
### Single Chr Visual ###
#########################

#Prepare Visual Objects
chrom = 3
chr.data <- data.table()
ids_to_vis <- which(colnames(adt.bychr) %in% high_qc_ids)
chr.data$TPM <- as.numeric(adt.bychr[c(chrom), ..ids_to_vis])
chr.data$VAR <- as.numeric(var_mat.bychr[c(chrom), ..ids_to_vis])
chr.data$DIFF <- as.numeric(varabs_mat.bychr[c(chrom), ..ids_to_vis])
chr.data$chr <- as.factor(rep(unique(adt.bychr$seqnames[c(chrom)]), length = length(chr.data$VAR)))
chr.data <- data.table(chr.data)

#Second Draft Plots
chr.visual <- ggplot(data = chr.data, mapping = aes(x = TPM, y = VAR, color = chr)) +
  geom_point(mapping = aes(x = TPM, y = VAR, color = chr)) +
  #xlim(c(0,2)) + ylim(c(0,2)) +
  labs(x = "TPM Ratio", y = "SNP Count Ratio", title = sprintf("Expression Contour Plot | chr%s ", chrom))

####################
### Cluster Plot ###
####################
#shorten data name
df <- chr.data[, c("TPM", "VAR", "DIFF")]
if (length(which(df$TPM > 2)) > 0) {df <- df[-which(df$TPM > 2),]}

#Visualize optimal number of clusters
wss <- fviz_nbclust(df, kmeans, method = "wss")
silhouette <- fviz_nbclust(df, kmeans, method = "silhouette")

#Extract optimal number of means
optimal_means_data <- fviz_nbclust(df, kmeans, method = "silhouette")$data
optimal_means <-as.numeric(optimal_means_data$clusters[which.max(optimal_means_data$y)])

#Visualize cluster map using optimal number of means
cluster.data <- kmeans(df, centers = optimal_means, nstart = 25)
cluster.visual <- fviz_cluster(cluster.data, data = df, geom = "point", choose.vars = c("TPM", "VAR"))

grid.arrange(arrangeGrob(chr.visual), arrangeGrob(cluster.visual, silhouette, ncol=1, nrow=2), widths=c(2,1))
