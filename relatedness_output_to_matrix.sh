#!/bin/bash -l

infile="1224samples_BiSNPs_MAC_DP_Fmissing_HETall_5_colonies_rel.relatedness"
outfile="1224samples_BiSNPs_MAC_DP_Fmissing_HETall_5_colonies_rel.relatedness_matrix"
>$outfile

printf "indv_vs_indv\t" > $outfile

cat $infile | grep -w "^P27152_381" | awk '{print $2}' | tr '\n' '\t' >> $outfile
printf "\n" >> $outfile

dummy=""

for indv in `cat $infile | tail -n +2 | awk '{print $1}' | uniq`; do
	rel_data=`cat $infile | grep -w "^"$indv | awk '{print $3}' | tr '\n' '\t'`
	printf "%s\t%s%s\n" "$indv" "$dummy" "$rel_data" >> $outfile
	dummy=`printf "%s0\t" "$dummy"`
done


