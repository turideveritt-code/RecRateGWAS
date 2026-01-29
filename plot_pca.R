rm(list=ls())

# ***********should exclude high missingness colonies?

setwd("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/filt_input_YAPP")

eigvec <- read.table("concat_chr_GQ20_Fmissing_hetALL.vcf.gz_plink_pca.eigenvec_for_plot",header=FALSE)
eigval <- read.table("concat_chr_GQ20_Fmissing_hetALL.vcf.gz_plink_pca.eigenval",header=FALSE)
PC1 <- eigval$V1[1]/sum(eigval$V1)
PC2 <- eigval$V1[2]/sum(eigval$V1)

col_vec <- c('grey57','cadetblue3','aquamarine4', 'darkgreen', 'olivedrab3', 'khaki','darkgoldenrod1', 'tomato', 'red3', 'violetred4','darkmagenta','deeppink2', 'plum3','deepskyblue1','slateblue','navy')

colony_list <- unique(eigvec$V1)

eigvec$color <- rep("gray50",nrow(eigvec))
eigvec$IDX <- 1:nrow(eigvec)

AA_list <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/MLH1_snp_GT_AA_2.txt",header=FALSE)
RR_list <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/MLH1_snp_GT_RR_2.txt",header=FALSE)

CO_list0 <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/filt_input_YAPP/YAPP_output/num_CO_sorted.txt",header=FALSE) # colony indv numCO
CO_list <- subset(CO_list0, V3>10 & V3 < 100)
meanCO <- mean(CO_list$V3,na.rm=TRUE)

CO_high <- subset(CO_list,V3>meanCO)
CO_low <- subset(CO_list,V3<=meanCO)

to_update <- subset(eigvec,V2 %in% CO_high$V2)
eigvec$color[to_update$IDX] <- "red3"

to_update <- subset(eigvec,V2 %in% CO_low$V2)
eigvec$color[to_update$IDX] <- "blue3"

# i <- 1
# for (colony in colony_list){
#   to_update <- subset(eigvec,eigvec$V1==colony)
#   eigvec$color[to_update$IDX] <- col_vec[i]
#   i <- i+1
#   if (i>length(col_vec)){i <- 1}
# }
plot(eigvec$V3, eigvec$V4, type="p", pch=20, cex=0.8, col=eigvec$color, xlab=paste("PC1 = ",100*round(PC1,digits=3),"%",sep=""), ylab=paste("PC2 = ",100*round(PC2,digits=3),"%",sep=""),main="Blue = low CO, Red = high CO")

write.table(CO_high$V2,file="drones_high_CO.txt",append=FALSE,quote=FALSE,row.names=FALSE, col.names=FALSE,sep="\n")
write.table(CO_low$V2,file="drones_low_CO.txt",append=FALSE,quote=FALSE,row.names=FALSE, col.names=FALSE,sep="\n")

eigvec$color <- "NA"
num_plots <- ceiling(length(colony_list)/length(col_vec))
par(mfrow=c(2,2))

ii <- 1
for (i in 1:num_plots){
  for (j in 1:length(col_vec)){
    if (ii <= length(colony_list)){
      to_update <- subset(eigvec,eigvec$V1==colony_list[ii])
      eigvec$color[to_update$IDX] <- col_vec[j]
      ii <- ii+1
    }
  }
  plot(eigvec$V3, eigvec$V4, type="p", pch=20, cex=0.8, col=eigvec$color, xlab=paste("PC1 = ",100*round(PC1,digits=3),"%",sep=""), ylab=paste("PC2 = ",100*round(PC2,digits=3),"%",sep=""))
  legend("topright",colony_list[((i-1)*length(col_vec)+1):(ii-1)],col=col_vec,pch=20,cex=0.6)
  eigvec$color <- "NA"
}

plot(1:nrow(eigval),eigval$V1/sum(eigval$V1),xlab="PC",ylab="Fraction explained")