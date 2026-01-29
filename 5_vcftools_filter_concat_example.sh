#!/bin/bash -l

#SBATCH -A snic2021-23-365
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J vcftools
#SBATCH --mail-type=all
#SBATCH --mail-user=example_email

# Main directories as variables
# Path to dir with all vcf-files
SEQDIR="/proj/snic2020-6-58/private/seq_data/drones_02974/P14101/contigs_all_no_136"
# Directory for output
OUTDIR=$SEQDIR/vcftools_output

cd $SEQDIR

module load bioinfo-tools
module load vcftools

for vcf in `ls *_output.vcf.hard.filtered.vcf`; do

	vcftools --vcf $vcf --max-missing 0.9 --min-alleles 2 --max-alleles 2 --remove-indels --remove-filtered-all --recode --out $OUTDIR"/"$vcf"_SNP_biallelic"
    
done

cd $OUTDIR

file_list=`ls *_output.vcf.hard.filtered.vcf_SNP_biallelic.recode.vcf`

module load bcftools

bcftools concat -O z -o $OUTDIR"/MR_SNP_biallelic_concat.vcf.gz" $file_list


