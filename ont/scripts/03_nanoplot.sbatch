#!/bin/bash
#SBATCH -J nanoplot
#SBATCH -o logs/nanoplot-%j.out
#SBATCH -c 4
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=00-01:00
#SBATCH --mem=5G
#SBATCH --array=1-5

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate nanoplot

# Define parameters
THREADS=8
INPUT_DIR="03_concatenate"
OUTPUT_DIR="04_nanoplot"

# Capture barcode name
BARCODE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f3)
echo $BARCODE

# Run NanoPlot
NanoPlot -t $THREADS \
	 --fastq $INPUT_DIR/"$BARCODE"*.fastq \
	 -o $OUTPUT_DIR/$BARCODE

