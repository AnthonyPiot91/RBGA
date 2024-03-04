#!/bin/bash
#SBATCH -J nanoplot
#SBATCH -o logs/nanoplot-%j.out
#SBATCH -c 4
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=00-01:00
#SBATCH --mem=5G
#SBATCH --array=1-8

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate nanoplot

# Define paramaters
THREADS=8
INPUT_DIR="06_chopper"
OUTPUT_DIR="07_nanoplot"
PARAM="2000_15"

# Capture strain name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

# Run NanoPlot
NanoPlot -t $THREADS \
	 --fastq $INPUT_DIR/"$STRAIN"_trimmed_"$PARAM".fastq \
	 -o $OUTPUT_DIR/$PARAM/$STRAIN/

