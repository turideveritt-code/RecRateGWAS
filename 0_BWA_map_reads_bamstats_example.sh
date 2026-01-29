#!/bin/bash -l

#SBATCH -A snic2021-22-195
#SBATCH -p core
#SBATCH -n 12
#SBATCH -t 30:00:00
#SBATCH -J bwa_mapping
#SBATCH --mail-type=all
#SBATCH --mail-user=example@email.com


# Main directories as variables
# SEQDIR assumed to contain all fastq-files

# Output: 	Sorted, indexed BAM-file for each sample
# 		BAM-stats for all samples written to the same file (bamtools_stats.txt in SEQDIR)

# Update this with the full path to the directory where you have the sample fastq files
SEQDIR="/proj/snic2020-6-58/private/seq_data/Am_drones_rerun_09-2020/P17406"

# Reference sequence
# Update this with the full path to the directory where you have the ref genome
REFDIR="/proj/snic2020-6-58/uppstore2017143/b2011048_nobackup/reference_genomes/Amel_HAv3.1"
# Update this with the name of the ref fasta file
REF=$REFDIR/GCF_003254395.2_Amel_HAv3.1_genomic.fna

map_reads () {

	# Load the bwa module to make bwa available
	# Load samtools, bamtools

	module load bioinfo-tools
	module load bwa
	module load samtools
	module load bamtools
	
	# Enter directory with data...
	
	cd $SEQDIR/
	
	# Create file for bamtools stats
	STAT_FILE=$SEQDIR/"bamtools_stats.txt"
	>$STAT_FILE
	
	for SAMPLE in `ls *_1.fastq | sed s/"_1.fastq"/""/g`; do #list all the sample names without the file extensions

		    			echo "MAPPING $SAMPLE ..."
                    
                    			# Run bwa across all allocated cores with the new mem
                    			# algorithm, and pipe the output to samtools to make a
                    			# sorted BAM file
                    
                    			time bwa mem -t 12 $REF \
						${SAMPLE}"_1.fastq" \
						${SAMPLE}"_2.fastq" \
						2> ${SAMPLE}.bam.bwa.log \
						| samtools view -b | samtools sort -o ${SAMPLE}.sorted.bam -O BAM
                    
                    			# Index the BAM file
                    
                    			samtools index ${SAMPLE}.sorted.bam

				# If sorted BAM-file was created in previous step, or existed from before: extract the bamstats and append to file
				if [ -e ${SAMPLE}.sorted.bam ]; then

					echo "Running bamtools stats for ${SAMPLE}.sorted.bam"
					echo "${SAMPLE}" >> $STAT_FILE
					bamtools stats -in  ${SAMPLE}.sorted.bam >> $STAT_FILE

				else

					echo "BAM-file ${SAMPLE}.sorted.bam not found. No stats produced."

				fi
        
	
	done

}

# Run the function map_reads above by uncommenting and execute the script

map_reads






