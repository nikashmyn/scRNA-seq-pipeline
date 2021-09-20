#################################################################
### Cluster/Contour Map Function (TPM Ratio vs Variant Ratio) ###
#################################################################

ClusterMap3D_byfamily <- function(adt.bychr, var_mat_sep.bychr, abs_allele_diff_mat.bychr, myids, myfamily, exclude_chrs) {
  
  ###################
  ### Contour Map ###
  ###################
  
  #Set seed for reproducability
  set.seed(123)
  
  #Prepare Visual Objects
  visual.data <- data.table(TPM = as.numeric(flatten(adt.bychr[,..myids])))
  visual.data$VAR <- as.numeric(flatten(var_mat_sep.bychr[,..myids]))
  visual.data$ABS <- as.numeric(flatten(abs_allele_diff_mat.bychr[,..myids]))
  visual.data$chr <- as.factor(rep(unique(adt.bychr$seqnames), length = length(visual.data$VAR)))
  visual.data$cell <- rep(myids, each=23-length(exclude_chrs))
  
  ####################
  ### Cluster Plot ###
  ####################
  #shorten data name
  df <- visual.data[, c("TPM", "VAR", "ABS")]
  
  #Visualize optimal number of clusters
  wss <- fviz_nbclust(df, kmeans, method = "wss")
  silhouette <- fviz_nbclust(df, kmeans, method = "silhouette")
  
  #Extract optimal number of means
  optimal_means_data <- fviz_nbclust(df, kmeans, method = "silhouette")$data
  optimal_means <-as.numeric(optimal_means_data$clusters[which.max(optimal_means_data$y)])
  
  #Visualize cluster map using optimal number of means
  cluster.data <- kmeans(df, centers = optimal_means, nstart = 25)
  #cluster.data <- kmeans(df, centers = 3, nstart = 25)
  
  ####################
  ### Final Visual ###
  ####################
  
  #TPM/VAR plane plot
  VAR.visual <- ggplot(data = visual.data, mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
    geom_point(mapping = aes(x = TPM, y = VAR, color = chr, shape = cell)) +
    xlim(c(0,2)) + ylim(c(0,2)) +
    labs(x = "TPM Ratio", y = "SNP Count Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily))
  
  #TPM/Diff plane plot
  DIFF.visual <- ggplot(data = visual.data, mapping = aes(x = TPM, y = ABS, color = chr, shape = cell)) +
    geom_point(mapping = aes(x = TPM, y = ABS, color = chr, shape = cell)) +
    theme(legend.position = "none") +
    #xlim(c(0,2)) + ylim(c(0,2)) +
    labs(x = "TPM Ratio", y = "Allelic Fraction Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily))
  
  #Cluster visual
  cluster.visual <- fviz_cluster(cluster.data, data = df, geom = "point", choose.vars = c("TPM", "VAR"))
  
  grid.arrange(arrangeGrob(VAR.visual), arrangeGrob(cluster.visual, DIFF.visual, ncol=1, nrow=2), widths=c(2,1))

}
  

ClusterMap3D_bychr <- function(D1, D2, D3, chr, high_qc_ids) {
  
  #########################
  ### Single Chr Visual ###
  #########################
  
  #Prepare Visual Objects
  chrom = 3
  chr.data <- data.table()
  ids_to_vis <- which(colnames(adt.bychr) %in% high_qc_ids)
  chr.data$TPM <- as.numeric(adt.bychr[c(chrom), ..ids_to_vis])
  chr.data$VAR <- as.numeric(var_mat.bychr[c(chrom), ..ids_to_vis])
  chr.data$ABS <- as.numeric(varabs_mat.bychr[c(chrom), ..ids_to_vis])
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
  df <- chr.data[, c("TPM", "VAR", "ABS")]
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
  
}

ClusterMap3D_Groups <- function(adt.bychr, var_mat_sep.bychr, abs_allele_diff_mat.bychr, myids, myfamily, mychrs) {
  
  #########################
  ### Single Chr Visual ###
  #########################
  
  #Prepare Visual Objects
  chr.data <- data.table() 
  chr.data$chr <- as.factor(as.numeric(mychrs))
  chr.data$TPM <- as.numeric(sapply(1:length(myids), function (x) { myid <- myids[x]; adt.bychr[,..myid][c(which(adt.bychr$seqnames == mychrs[x]))] }))
  chr.data$VAR <- as.numeric(sapply(1:length(myids), function (x) { myid <- myids[x]; var_mat_sep.bychr[,..myid][c(which(adt.bychr$seqnames == mychrs[x]))] }))
  chr.data$ABS <- as.numeric(sapply(1:length(myids), function (x) { myid <- myids[x]; abs_allele_diff_mat.bychr[,..myid][c(which(adt.bychr$seqnames == mychrs[x]))] }))
  chr.data <- data.table(chr.data)
  
  print(chr.data)
  
  #Second Draft Plots
  chr.visual <- ggplot(data = chr.data, mapping = aes(x = TPM, y = VAR, color = chr)) +
    geom_point(mapping = aes(x = TPM, y = VAR, color = chr)) +
    #xlim(c(0,2)) + ylim(c(0,2)) +
    labs(x = "TPM Ratio", y = "SNP Count Ratio", title = sprintf("Expression Contour Plot | %s ", myfamily))
  
  ####################
  ### Cluster Plot ###
  ####################
  #shorten data name
  df <- chr.data[, c("TPM", "VAR", "ABS")]
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
  
}