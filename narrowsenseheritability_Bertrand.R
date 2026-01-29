rm(list=ls())

library(MM4LMM)
library(dplyr)

setwd("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/final_metadata")

# import r_intra data, already filtered for quality and recombination outliers

r_intra_init0 <- read.table("/Users/turev988/Sync/HoneyBees/drones_GWAS/final_dataset/r_intra_CO_interference/r_intra_per_drone_filt.txt",header=TRUE) #parent  offspring       r_intra colony_id       total_CO

### ******************* added 2026-01-07, very minor changes to results so not updated in manuscript
replicates_exclude <- c("VSH_004_8v2","VSH_004_7v2","VSH_004_6v2","VSH_004_5v2","VSH_004_4v2","VSH_004_3v2","VSH_004_2v2","VSH_004_1v2","PK_02II_8v2","PK_02II_7v2","PK_02II_5v2","PK_02II_4v2","PK_02II_3v2","PK_02II_2v2","PK_02II_1v2","IA_OS24_10v2","HJ_H3_15v2")
r_intra_init <- subset(r_intra_init0,!(offspring %in% replicates_exclude))
### ****************************************

# import observed recomb data
# not filtered, but filtering by combining with r_intra
obs_recomb00 <- read.table("all_chr_all_drones_CO_per_bp_corrected_genome_size2.txt",header=TRUE) #colony	indv	num_CO	indv_same	corrected_genome_bp	CO_per_bp

r_intra <- left_join(r_intra_init,obs_recomb00,by=join_by(offspring==indv)) # contains both r_intra and recomb rate

# To control for CO-freq when analysing r_intra:
# linear regression of r_intra on CO-freq, then take the residuals of that as the r_intra phenotype controlled for CO-freq. 
lm_intra_recomb <- lm(r_intra ~ CO_per_bp, data=r_intra)
r_intra$linreg_residuals <- lm_intra_recomb$residuals

num_drones <- nrow(r_intra)
queen_list <- unique(r_intra$colony_id)
num_queens <- length(queen_list)

##################

# import GRM
#grm_table <- read.table("Demetris_matrix.cXX.txt",header=FALSE)
#grm_table <- read.table("LDAKgmat.grm.raw",header=FALSE)
#grm_queens_all <- read.table("Demetris_queens.fam_original",header=FALSE) #colony names corresponding to GRM matrix

grm_table <- read.table("unimputed_data/LDAKgmat.grm.raw",header=FALSE)
grm_queens_all <- read.table("unimputed_data/all_queen_phased.fam",header=FALSE) #colony names corresponding to GRM matrix

rownames(grm_table) <- grm_queens_all$V1 #rownames and colnames for relatedness plot
colnames(grm_table) <- grm_queens_all$V1
grm_table$queens <- grm_queens_all$V1
grm_table$IDX <- 1:nrow(grm_table)

# only include the same colonies in the GRM as in the (already filtered) r_intra table
grm_include <- subset(grm_table,queens %in% queen_list)

grm_table_filt <- grm_table[grm_include$IDX,grm_include$IDX]
grm_matrix <- as.matrix(grm_table_filt)
#heatmap(grm_matrix,Rowv=NA,Colv=NA,scale="none",main="LDAK")

#sort queens in r_intra in same order as in grm_include$queens

r_intra$dummy <- rep(0,nrow(r_intra))
r_intra$IDX <- 1:nrow(r_intra)
for (i in 1:nrow(grm_include)){
  queen_i <- grm_include$queens[i]
  ri_subset <- subset(r_intra,colony_id==queen_i)
  r_intra$dummy[ri_subset$IDX] <- i
}

r_intra_order <- order(r_intra$dummy)

r_intra_sorted <- r_intra[r_intra_order,]
queen_list_sorted <- unique(r_intra_sorted$colony_id)

r_intra_sorted$IDX <- 1:nrow(r_intra_sorted)

# create all input matrices for MMEst

VL_m <- diag(num_queens)
VL_a <- grm_matrix #diag(num_queens) #GRM - to be added
VL_e <- diag(num_drones)

ZL_m <- diag(0,num_drones,num_queens)
ZL_e <- diag(num_drones)

for (i in 1:num_queens){
  q_i <- queen_list_sorted[i]
  q_data <- subset(r_intra_sorted, colony_id==q_i)
  ZL_m[q_data$IDX,i] <- 1
}

## Just a trick to create these design matrices
## make sure to order the factor level as you want
r_intra_sorted$colony_id = factor(r_intra_sorted$colony_id,levels=queen_list_sorted)
## use model.matrix to build the design matrix with correct contrast for a random effect
ZL_m_mm = model.matrix( ~ 0 + r_intra_sorted$colony_id)
## they are the same matrices
sum(ZL_m_mm!=ZL_m) # 0

ZL_m = ZL_m_mm
ZL_a <- ZL_m

VL <- list(VL_m,VL_a,VL_e)
ZL <- list(ZL_m,ZL_a,ZL_e)

# Just scale the phenotype to a Gaussian distribution to avoid numerical problems
r_intra_sorted$Y_COfreq = qqnorm(r_intra_sorted$CO_per_bp)$x 
r_intra_sorted$Y_r_intra = qqnorm(r_intra_sorted$r_intra)$x 
r_intra_sorted$Y_r_intra_resid = qqnorm(r_intra_sorted$linreg_residuals)$x 

# Animal Model fitted with ML4LMM
reml_COfreq <- MMEst(Y=r_intra_sorted$Y_COfreq,VarList=VL,ZList=ZL,Method="Reml")
reml_r_intra <- MMEst(Y=r_intra_sorted$Y_r_intra,VarList=VL,ZList=ZL,Method="Reml")
reml_r_intra_resid <- MMEst(Y=r_intra_sorted$Y_r_intra_resid,VarList=VL,ZList=ZL,Method="Reml")

###

sigma_m <- as.numeric(reml_COfreq$NullModel$Sigma2[1])
sigma_a <- as.numeric(reml_COfreq$NullModel$Sigma2[2])
sigma_e <- as.numeric(reml_COfreq$NullModel$Sigma2[3])

H2_broad_COfreq <- (sigma_m+sigma_a)/(sigma_m+sigma_a+sigma_e)
h2_narrow_COfreq <- (sigma_a)/(sigma_m+sigma_a+sigma_e)
perm_env_effect_COfreq = sigma_m/(sigma_m+sigma_a+sigma_e)

###

sigma_m <- as.numeric(reml_r_intra$NullModel$Sigma2[1])
sigma_a <- as.numeric(reml_r_intra$NullModel$Sigma2[2])
sigma_e <- as.numeric(reml_r_intra$NullModel$Sigma2[3])

H2_broad_r_intra <- (sigma_m+sigma_a)/(sigma_m+sigma_a+sigma_e)
h2_narrow_r_intra <- (sigma_a)/(sigma_m+sigma_a+sigma_e)
perm_env_effect_r_intra = sigma_m/(sigma_m+sigma_a+sigma_e)

###

sigma_m <- as.numeric(reml_r_intra_resid$NullModel$Sigma2[1])
sigma_a <- as.numeric(reml_r_intra_resid$NullModel$Sigma2[2])
sigma_e <- as.numeric(reml_r_intra_resid$NullModel$Sigma2[3])

H2_broad_r_intra_resid <- (sigma_m+sigma_a)/(sigma_m+sigma_a+sigma_e)
h2_narrow_r_intra_resid <- (sigma_a)/(sigma_m+sigma_a+sigma_e)
perm_env_effect_r_intra_resid = sigma_m/(sigma_m+sigma_a+sigma_e)


####################### main analysis done, extra stuff below
## Queen effect, anova
mod.fix = lm( Y ~colony, data = obs_recomb_sorted)
anova(mod.fix)
## from this we can estimate the variance components of a queen model
## I copy pasted the values from the anova table
SSE = 718.77
SSA = 788.93
m = 184
n = 1325

Ve_aov = SSE/n
Vq_aov = (SSA - (m-1)*Ve_aov)/n

H2_broad_aov = Vq_aov/(Ve_aov+Vq_aov)
H2_broad_aov

## Letś try with nlme
library(nlme)
lme.data = groupedData(Y ~1|colony,obs_recomb_sorted) # |> select(Y,colony)
mod.ran =  lme(Y ~ 1, data = lme.data, random = ~1 | colony )
VC = as.matrix(VarCorr(mod.ran))
Vq = as.double(VC[1,1])
Ve = as.double(VC[2,1])
## This is "Broad Sense Heritability" = random queen effect model
H2_broad_reml = Vq/(Vq+Ve)
H2_broad_reml
