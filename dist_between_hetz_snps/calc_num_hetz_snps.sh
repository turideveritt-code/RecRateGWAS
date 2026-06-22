#!/bin/bash -l

# make one separate vcf-file for each queen, which only includes SNPs that are heterozygous for that queen
# write a list of all queens and their numbers of heterozygous snps

vcf="concat_queen_biallelic_rm_bad_indv_mac1_RIC.vcf"

>hetz_snps_per_queen.txt

for queen in `cat "queen_samp_names.txt"`; do
	vcftools --vcf $vcf --indv $queen --mac 1 --recode --out $vcf"_"$queen
	hetz_count=`cat $vcf"_"$queen".recode.vcf" | grep -v "^#" | wc -l`
	printf "%s\t%d\n" $queen $hetz_count >> hetz_snps_per_queen.txt
done


