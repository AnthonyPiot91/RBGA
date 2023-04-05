#!/bin/bash
#SBATCH -J minimap2
#SBATCH -o logs/minimap2-%j.out
#SBATCH -c 6
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=3G
#SBATCH --array=1-5

# Load conda environment
source /mnt/ibis/rclevesq/software/miniconda3/anpio1/etc/profile.d/conda.sh
conda activate minimap2

# Define paramaters
THREADS=12
INPUT_DIR_FASTQ="minimap2/chromosome"
OUTPUT_DIR="minimap2/contigs_wo_chrom"

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
REF="canu_plasmids/$STRAIN/$STRAIN.contigs.fasta"
echo $REF

minimap2 --version

# Run minimap2
minimap2 -a \
	 -x map-ont \
	 $REF \
	 $INPUT_DIR_FASTQ/"$STRAIN"_chromosome_unmapped.fastq \
	 -o $OUTPUT_DIR/"$STRAIN"_contigs_wo_chrom.sam

# Deactivate minimap2 environment to use samtools
conda deactivate

samtools --version | sed -n '1,3p'

# Extract unmapped reads with samtools
samtools view -bf 4 \
		$OUTPUT_DIR/"$STRAIN"_contigs_wo_chrom.sam \
		> $OUTPUT_DIR/"$STRAIN"_contigs_wo_chrom_unmapped.bam

# Convert bam file to fastq
samtools fastq $OUTPUT_DIR/"$STRAIN"_contigs_wo_chrom_unmapped.bam \
	     > $OUTPUT_DIR/"$STRAIN"_contigs_wo_chrom_unmapped.fastq
