#!/bin/bash
#SBATCH -J flye
#SBATCH -o logs/flye-%j.out
#SBATCH -c 6
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=00-01:30
#SBATCH --mem=10G
#SBATCH --array=1-6

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate flye

# Define paramaters
THREADS=8
GENOME_SIZE=1500000
INPUT_DIR="06_chopper"
OUTPUT_DIR="10_flye"

# Capture strain and fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
INPUT_FILE_PATH="$INPUT_DIR/"$STRAIN"*.fastq"
echo $INPUT_FILE_PATH

# Run flye
flye 	--nano-raw \
	$INPUT_FILE_PATH \
	--genome-size $GENOME_SIZE \
	--threads $THREADS \
	--out-dir $OUTPUT_DIR/"$STRAIN"
