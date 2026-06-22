rm(list=ls())

# sum of hetz SNPs per queen
snps_per_queen <- read.table("hetz_snps_per_queen.txt",header=FALSE)
par(mfrow=c(1,1))
barplot(snps_per_queen$V2,names.arg=snps_per_queen$V1,las=2,cex.names=0.7)
genome_size <- 225200000
mean_dist_betw_snps <- (nrow(snps_per_queen)*genome_size)/sum(snps_per_queen$V2)

# list of all distances between consecutive hetz SNPs (per queen and chrom)
all_dists <- read.table("snp_dists_for_hist.txt",header=FALSE)
par(mfrow=c(1,3))
hist(all_dists$V1, main="Dists between hetz SNPs")
hist(all_dists$V1,xlim=c(0,4e+03),breaks=80000, main="Dists between hetz SNPs")
summary(all_dists$V1)

# compare to distances between consecutive COs (per drone and chrom)

CO_data <- read.table("concat_ID_queen_biallelic_yapp_recombinations_filtered.txt",header=TRUE)
CO_data$mid_pos <- (CO_data$left + CO_data$right)/2

CO_dists_all <- rep(0,nrow(CO_data))

chr_list <- unique(CO_data$chrom)
drone_list <- unique(CO_data$offspring)


i <- 1
for (drone in drone_list){
  for (chrom_i in chr_list){
    CO_subset <- subset(CO_data,chrom==chrom_i & offspring==drone)
    #print(nrow(CO_subset))
    if (nrow(CO_subset)>0){
      CO_dists_all[i] <- CO_subset$mid_pos[1]
      i <- i+1
      if (nrow(CO_subset)>1){
        for (j in 2:nrow(CO_subset)){
          CO_dists_all[i] <- CO_subset$mid_pos[j] - CO_subset$mid_pos[(j-1)]
          i <- i+1
        }
      }
    }
  }
}

hist(CO_dists_all,main="Dists between consecutive COs")
summary(CO_dists_all)

#heterozygous SNPs
summary(all_dists$V1)
