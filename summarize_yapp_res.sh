#!/bin/bash -l

YAPPfile="merged_yapp_recombinations.txt"

out_file=$YAPPfile"_matrix"

>$out_file

chr_list=("NC_037638.1" "NC_037639.1" "NC_037640.1" "NC_037641.1" "NC_037642.1" "NC_037643.1" "NC_037644.1" "NC_037645.1" "NC_037646.1" "NC_037647.1" "NC_037648.1" "NC_037649.1" "NC_037650.1" "NC_037651.1" "NC_037652.1" "NC_037653.1")

printf "colony\tindv" > $out_file
for chr in `echo ${chr_list[*]}`; do
	printf "\t%s" $chr >> $out_file
done
printf "\tsumCO\n" >> $out_file

indv_list=`cat $YAPPfile | awk '{print $3}' | grep -v "Agroup" | grep -v "Cgroup" | grep -v "Mgroup" | grep -v "offspring" | sort | uniq`

for indv in `echo ${indv_list[*]}`; do
	colony=`cat $YAPPfile | grep -w $indv | head -n 1 | awk '{ print $1 }' | sed s/"_queen"/""/g`
	printf "%s\t%s" $colony $indv >> $out_file
	for chr in `echo ${chr_list[*]}`; do
		num_CO=`cat $YAPPfile | grep -w $indv | grep -w $chr | wc -l`
		printf "\t%d" $num_CO >> $out_file
	done
	tot_CO=`cat $YAPPfile | grep -w $indv | wc -l`
	printf "\t%d\n" $tot_CO >> $out_file
done


