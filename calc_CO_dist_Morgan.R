rm(list=ls())

setwd("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/filt_input_YAPP/YAPP_output")

CO_all_concat <- read.table("concat_ID_queen_biallelic_yapp_recombinations_filtered.txt",header=TRUE)

### look at distribution of interval distances
CO_all_concat$interv_dist <- CO_all_concat$right - CO_all_concat$left
hist(CO_all_concat$interv_dist)
dist_mean <- mean(CO_all_concat$interv_dist,na.rm=TRUE)
dist_sd <- sd(CO_all_concat$interv_dist,na.rm=TRUE)

###

CO_all_concat$midpoint <- (CO_all_concat$left + CO_all_concat$right)/2

CO_pos_order <- order(CO_all_concat$chrom,CO_all_concat$midpoint)

CO_all_sorted_pos <- CO_all_concat[CO_pos_order,]

CO_all_concat$cM_to_previous_CO <- rep(0,nrow(CO_all_concat))

num_drones <- length(unique(sort(CO_all_concat$offspring)))

for (i in 1:nrow(CO_all_concat)){
  chr_i <- CO_all_concat$chrom[i]
  indv_i <- CO_all_concat$offspring[i]
  bp_i <- CO_all_concat$midpoint[i]
  if (i>1){
    if (chr_i == CO_all_concat$chrom[(i-1)] & indv_i == CO_all_concat$offspring[(i-1)]){
      bp_before <- CO_all_concat$midpoint[(i-1)]
      #print(bp_before)
    }
    else {
      bp_before <- 0
      #print(bp_before)
    }
  }
  else {
    bp_before <- 0
  }
  CO_all_concat$cM_to_previous_CO[i] <- nrow(subset(CO_all_sorted_pos, CO_all_sorted_pos$chrom==chr_i & CO_all_sorted_pos$midpoint>bp_before & CO_all_sorted_pos$midpoint<bp_i))*100/num_drones
}

#write.table(CO_all_concat,file="/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/r_intra_CO_interference/concat_ID_queen_biallelic_yapp_recombinations_filtered_CO_dist_cM.txt",col.names=TRUE,row.names = FALSE,quote=FALSE,sep="\t",append=FALSE)