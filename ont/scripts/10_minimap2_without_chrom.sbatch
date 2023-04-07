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
INPUT_DIR_FASTQ="10_minimap2/chromosome"
OUTPUT_DIR="10_minimap2/without_chrom"
REF_DIR="09_canu/without_chrom"

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
REF="$REF_DIR/$STRAIN/$STRAIN.contigs.fasta"
echo $REF

minimap2 --version

# Run minimap2
minimap2 -a \
	 -x map-ont \
	 $REF \
	 $INPUT_DIR_FASTQ/"$STRAIN"_chromosome_unmapped.fastq \
	 -o $OUTPUT_DIR/"$STRAIN"_contigs_wo_chrom.sam
