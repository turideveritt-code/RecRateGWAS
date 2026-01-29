rm(list=ls())
setwd("/Users/turev988/Sync/HoneyBees/drones_GWAS/troubleshooting/test_colonies/replace_missing_by_dots/yapp_v2")

rel_data <- read.table("1224samples_BiSNPs_MAC_DP_Fmissing_HETall_5_colonies_rel.relatedness_matrix",header=TRUE, row.names=1)

heatmap(as.matrix(rel_data),Rowv=NA,Colv=NA,scale="none",col=hcl.colors(20, palette = "reds 2", alpha = 1,rev=TRUE))

legend(x="topleft", legend=c("min", "x", "x", "x", "max"), fill=hcl.colors(5, palette = "reds 2", alpha = 1,rev=TRUE))
