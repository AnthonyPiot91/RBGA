#!/bin/bash
#SBATCH -J biopython
#SBATCH -o logs/biopython-%j.out
#SBATCH -c 1
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=3G
#SBATCH --array=1-5

# Load conda environment
source /mnt/ibis/rclevesq/software/miniconda3/anpio1/etc/profile.d/conda.sh

# Define paramaters
INPUT_DIR="09_canu/with_chrom"
OUTPUT_DIR="10_minimap2/chromosome"

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

# Run python script
python  scripts/sub_scripts/extract_longest_fasta_sequence.py \
	$INPUT_DIR/$STRAIN/$STRAIN.contigs.fasta \
	$OUTPUT_DIR/"$STRAIN"_chromosome.fasta
