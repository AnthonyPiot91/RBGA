#!/bin/bash
#SBATCH -J nanoQC
#SBATCH -o logs/nanoQC-%j.out
#SBATCH -c 1
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:20
#SBATCH --mem=6G
#SBATCH --array=1-8

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate RBGA_1

# Define paramaters
INPUT_DIR="06_chopper"
OUTPUT_DIR="08_nanoQC"
PARAM="2000_15"

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

# Create output subdir
if [ ! -d $OUTPUT_DIR/$PARAM/$STRAIN ]
	then mkdir -p $OUTPUT_DIR/$PARAM/$STRAIN
fi

# Run nanoQC
nanoQC -o $OUTPUT_DIR/$PARAM/$STRAIN $INPUT_DIR/"$STRAIN"_trimmed_"$PARAM".fastq
