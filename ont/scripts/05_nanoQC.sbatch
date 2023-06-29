#!/bin/bash
#SBATCH -J nanoQC
#SBATCH -o logs/nanoQC-%j.out
#SBATCH -c 1
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:20
#SBATCH --mem=4G
#SBATCH --array=1-5

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate RBGA_1

# Define paramaters
INPUT_DIR="03_concatenate"
OUTPUT_DIR="05_nanoQC"

# Capture fastq file name
BARCODE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f3)
echo $BARCODE

# Create output subdir
if [ ! -d $OUTPUT_DIR/$BARCODE ]
	then mkdir -p $OUTPUT_DIR/$BARCODE
fi

# Run nanoQC
nanoQC -o $OUTPUT_DIR/$BARCODE $INPUT_DIR/"$BARCODE"*.fastq
