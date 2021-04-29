#import gene list
full_genes <- read.csv("/projectnb/bf528/users/lava_lamp/project_4/marker_genes.csv")

for (i in 1:11) {
  df <- data.frame(full_genes$gene[full_genes$cluster == i])
  write.table(df, file = paste("/projectnb/bf528/users/lava_lamp/project_4/alec_pr4/cluster_",i,".txt",sep = ""),
              quote = F, row.names = F, col.names = F)
  assign(paste("cluster_",i,sep=""),df)
}

filtered_genes <- full_genes[full_genes$p_val_adj <= 0.05 & full_genes$avg_log2FC >= 1.5,]

for (i in 1:11) {
  df <- data.frame(filtered_genes$gene[filtered_genes$cluster == i])
  write.table(df, file = paste("/projectnb/bf528/users/lava_lamp/project_4/alec_pr4/filt_cluster_",i,".txt",sep = ""),
              quote = F, row.names = F, col.names = F)
  assign(paste("filtered_cluster_",i,sep=""),df)
}
