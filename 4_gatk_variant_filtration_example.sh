#!/bin/bash -l

#SBATCH -A snic2021-22-195
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J gatk_hard_filt
#SBATCH --mail-type=all
#SBATCH --mail-user=example@example.se

# Main directories as variables

# Path to where you have the vcf-files
SEQDIR="/proj/snic2020-6-58/private/seq_data/drones_02974/P14101/contigs_all_no_136"

# Reference sequence
# Path and file name for ref genome
REFDIR="/proj/snic2020-6-58/uppstore2017143/b2011048_nobackup/reference_genomes/Amel_HAv3.1"
REF=$REFDIR/GCF_003254395.2_Amel_HAv3.1_genomic.fna

function filter_variants() {

	module load bioinfo-tools
	module load GATK

	cd $SEQDIR

	for vcf in `ls *_output.vcf`; do

		time java -jar -Xmx4g -jar /sw/apps/bioinfo/GATK/3.7/GenomeAnalysisTK.jar \
			-T VariantFiltration \
			-R $REF \
			-V $vcf \
			--filterExpression "QD<2.0" \
	    		--filterName "QD" \
	    		--filterExpression "FS>60.0" \
	    		--filterName "FS" \
	    		--filterExpression "MQ<40.0" \
	    		--filterName "MQ" \
	    		--filterExpression "MQRankSum<-12.5" \
	    		--filterName "MQRankSum" \
	    		--filterExpression "ReadPosRankSum<-8.0" \
	    		--filterName "ReadPosRankSum" \
			--filterExpression "SOR>3.0" \
			--filterName "SOR" \
	    		-o ${vcf}.hard.filtered.vcf
    
	done

}

filter_variants











