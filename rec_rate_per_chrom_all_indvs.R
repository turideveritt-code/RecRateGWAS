rm(list=ls())

# Plot recomb rate per pos for all chromosomes
# Compare with LDAK weights - correlation?

num_indvs <- 1508

recomb_map <- read.table("/Users/turev988/Downloads/Supplementary_Table_4_recomb_per_SNP_all_concat.txt",header=TRUE)
chrom_list <- unique(recomb_map$chrom)

#weights_ldak <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/final_metadata/LDAKgmat_details_with_chrpos_sorted.txt",header=TRUE)
# weights_ldak <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/final_metadata/unimputed_data/LDAKgmat.grm.details_chr_pos",header=TRUE)
# weights_ldak$IDX <- 1:nrow(weights_ldak)
# weights_ldak$chr_name <- rep("NA",nrow(weights_ldak))

centrom <- read.table("/Users/turev988/Downloads/evad157_supplementary_data/supp_tables_S1_S6_postrev.csv",header=TRUE,sep=";")

par(mfrow=c(4,4))

j=1
for (chrom_i in chrom_list) {
  # chrom_weights <- subset(weights_ldak,Chromosome==chrom_i)
  chrom_data <- subset(recomb_map,chrom==chrom_i); 
  # weights_ldak$chr_name[chrom_weights$IDX] <- chrom_i 
  plot((chrom_data$start_pos+chrom_data$end_pos)/2,(chrom_data$CO_per_bp)*(10^8)/num_indvs,type="p",pch=20,cex=0.3,main=paste("Chr",j),xlab="pos",ylab="cM/Mb"); #,ylim=c(0,150));
  abline(v=centrom$start_pos[j],col="red3",lty=2);
  abline(v=centrom$end_pos[j],col="green3",lty=2);
  # plot(chrom_weights$Position,chrom_weights$Weight,pch=20,cex=0.3,col="gray40",main=chrom_i,ylab="Weights")
  # abline(v=centrom$start_pos[j],col="red3",lty=2);
  # abline(v=centrom$end_pos[j],col="green3",lty=2);
  j=j+1
}

# recomb_map$chr_pos <- paste(recomb_map$chrom,recomb_map$end_pos,sep="_")
# weights_ldak$chr_pos <- paste(weights_ldak$chr_name,weights_ldak$Position,sep="_")
# 
# recomb_weights_merge <- merge(recomb_map,weights_ldak,by.x="chr_pos",by.y="chr_pos")
# recomb_weights_merge$cMMb <- recomb_weights_merge$CO_per_bp*(10^8)/num_indvs
# 
# par(mfrow=c(1,1))
# 
# plot(recomb_weights_merge$Weight,recomb_weights_merge$cMMb,type="p",pch=20,cex=0.4)
# cor.test(recomb_weights_merge$Weight,recomb_weights_merge$cMMb,method="spearman",two.sided=TRUE)

