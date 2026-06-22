#!/bin/bash -l

# for each queen and each chromsome separately, calculate the distance between each pair of consecutive SNPs
# use the vcf-files which only include one queen and the heterozygous SNPs in that queen
# write all the distances to the same file

chr_list="chr_list.txt"
out_list="snp_dists_for_hist.txt"

>$out_list

for vcf in *_queen.recode.vcf; do
	for chrom in `cat $chr_list`; do
		cat $vcf | grep -v "^#" | grep $chrom | awk 'BEGIN {i_0=0} {print $2-i_0; i_0=$2}' >> $out_list
	done
done

