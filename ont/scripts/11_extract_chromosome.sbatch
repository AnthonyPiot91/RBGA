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
INPUT_DIR="canu"
OUTPUT_DIR="minimap2"

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
CHROM_REF="$OUTPUT_DIR/"$STRAIN"_chromosome.fasta"
echo $CHROM_REF

# Run python script
python  scripts/sub_scripts/extract_longest_fasta_sequence.py \
	canu/clean/$STRAIN/$STRAIN.contigs.fasta \
	minimap2/"$STRAIN"_chromosome.fasta
