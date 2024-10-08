#!/bin/bash
#SBATCH -J chopper
#SBATCH -o logs/chopper-%j.out
#SBATCH -c 8
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:30
#SBATCH --mem=1G
#SBATCH --array=2-5

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate chopper

# Define parameters
# For canu, MIN_QUAL and MIN_LEN can be kept at zero because this assembler will
# use in priority the longest reads as well as additional smaller reads to represent
# the whole genome if needed. Reads will then be corrected by canu.
# HEADCROP and TAILCROP should be set according the the 05_nanoQC results.
THREADS=8
MIN_QUAL=10
MIN_LEN=2000
HEADCROP=25
TAILCROP=25
INPUT_DIR="03_concatenate"
OUTPUT_DIR="06_chopper"

# Capture fastq file name
BARCODE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f3)
echo $BARCODE
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

# Run chopper
cat $INPUT_DIR/"$BARCODE".merged.fastq | \
chopper -t $THREADS \
	-q $MIN_QUAL \
	-l $MIN_LEN \
	--headcrop $HEADCROP \
	--tailcrop $TAILCROP \
> $OUTPUT_DIR/"$STRAIN"_trimmed_"$MIN_LEN"_"$MIN_QUAL".fastq
