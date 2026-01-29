#!/bin/bash -l

#SBATCH -A snic2021-22-195
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 50:00:00
#SBATCH -J haplotype_caller
#SBATCH --mail-type=all
#SBATCH --mail-user=example@example.se

# Main directories as variables

# Update to location of your bam-files
SEQDIR=/proj/uppstore2017143/b2011048_nobackup/seq_data/drones_02974/P14101

# This file is created to be used in the next step
SAMP_NAME_MAP=$SEQDIR/samp_name_map.txt
> $SAMP_NAME_MAP

# Update with the path to your ref genome
REFDIR=/proj/uppstore2017143/b2011048_nobackup/reference_genomes/Amel_HAv3.1
REF=$REFDIR/GCF_003254395.2_Amel_HAv3.1_genomic.fna

# Haploid or diploid genotyping?
PLOIDY=2


# This function is to run the jobs in parallel
function pwait() {
  while [ $(jobs -p | wc -l) -ge $1 ]; do
    sleep $2
  done
}

function call_variants() {

	module load bioinfo-tools
	module load GATK/4.0.8.0 #Need to use GATK4 for GenomicsDBImport

	cd $SEQDIR

	for sample in `ls MR_*.sorted.marked_dupes.rg.bam | sed s/".sorted.marked_dupes.rg.bam"/""/g`; do
		# The jobs are submitted in the background (with the "&" option below)
		# The pwait function prevents jobs from being submitted if the max number of jobs (20) are already running
		# Then it waits for 20 s and then checks again
		pwait 20 20s

		echo "Running HaplotypeCaller for $sample"

		gatk HaplotypeCaller \
		-R $REF \
		-I $sample".sorted.marked_dupes.rg.bam" \
		-O $sample".sorted.marked_dupes.rg.g.vcf.gz" \
		-ERC GVCF \
		-ploidy $PLOIDY &

		# list the sample names and directories in the samp_name_map - file
		printf "$sample\t`pwd`/$sample.sorted.marked_dupes.rg.g.vcf.gz\n" >> $SAMP_NAME_MAP
    
	done
	wait # wait until all jobs have finished

}

call_variants



