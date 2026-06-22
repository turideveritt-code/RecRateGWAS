rm(list=ls())

library(gaston)
library(stringr)


recomb_data00 <- read.table("all_chr_all_drones_CO_per_bp_corrected_genome_size2.txt",header=TRUE) #colony	indv	num_CO	indv_same	corrected_genome_bp	CO_per_bp
to_exclude <- read.table("excluded_samp.txt",sep="\t",header=TRUE) #colony  indv    nig_id  incl_excl       reason_excl
recomb_data0 <- subset(recomb_data00,!(indv %in% to_exclude$indv))
recomb_data <- subset(recomb_data0,num_CO<100) # remove outliers
num_drones <- nrow(recomb_data)

##########**********
## GRM from Gaston, include drones according to B's explanation

bed_test <- read.bed.matrix("concat_queen_biallelic_rm_bad_indv_mac1_RIC.vcf_queens_filt_chrom_int_plink")
standardize(bed_test) <- "p"
GRM_test <- GRM(bed_test)

queens_list <- bed_test@ped$id
num_queens <- length(queens_list)
# remove "_queen" from sample names
for (i in 1:num_queens){
  num_char <- str_length(queens_list[i])
  queens_list[i] <- substr(queens_list[i],1,(num_char-6))
}

Z <- matrix(data=0,nrow=num_queens,ncol=num_drones)
Y <- rep(-999,num_drones) # put phenotypes in same order as queens in GRM

# reorder other data as well to match the order in the GRM
recomb_data_new <- data.frame("colony"=rep("NA",num_drones),"indv"=rep("NA",num_drones),"num_CO"=rep(0,num_drones),"corrected_genome_bp"=rep(0,num_drones),"CO_per_bp"=rep(0,num_drones))

list_count <- 1
for (i in 1:num_queens){
  qi <- queens_list[i]
  colony_subset <- subset(recomb_data,colony==qi)
  num_in_qi <- nrow(colony_subset)
  Y[list_count:(list_count+num_in_qi-1)] <- colony_subset$CO_per_bp
  Z[i,list_count:(list_count+num_in_qi-1)] <- 1
  recomb_data_new[list_count:(list_count+num_in_qi-1),] <- colony_subset[,c("colony","indv","num_CO","corrected_genome_bp","CO_per_bp")]
  list_count <- list_count + num_in_qi
}

GRM_drones <- t(Z) %*% GRM_test %*% Z
#heatmap(GRM_test,Rowv=NA,Colv=NA,scale="none",main="GRM queens")
#heatmap(GRM_drones,Rowv=NA,Colv=NA,scale="none",main="GRM drones")

# scale phenotypes (Y) by 10^6 in order to avoid numerical issues
gaston_output_drones <- lmm.aireml(Y*(10^6),K=GRM_drones,verbose=TRUE)

plot(gaston_output_drones$BLUP_omega,Y,type="p",pch=20,cex=0.7,main="GRM from Gaston, t(Z)*GRM*Z")

recomb_data_new$blup_gaston <- gaston_output_drones$BLUP_omega

#write.table(recomb_data_new,file="drones_GWAS/final_dataset/redoing_gwas/blups_gaston_updatedApril26/blups_gaston_per_drone.tsv",row.names=FALSE,col.names = TRUE,quote = FALSE,sep="\t",append=FALSE)

# The drones have the same BLUP value within each colony
blup_per_queen <- rep(0,num_queens)

for (i in 1:num_queens){
  qi <- queens_list[i]
  colony_subset <- subset(recomb_data_new,colony==qi)
  blup_per_queen[i] <- mean(colony_subset$blup_gaston,na.rm=TRUE)
}

queen_blup_mean_vals <- data.frame("colony"=queens_list,"BLUP_mean"=blup_per_queen)

#write.table(queen_blup_mean_vals,file="drones_GWAS/final_dataset/redoing_gwas/blups_gaston_updatedApril26/blups_gaston_per_colony.tsv",row.names=FALSE,col.names = TRUE,quote = FALSE,sep="\t",append=FALSE)

