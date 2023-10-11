#!/bin/bash
#SBATCH -J fastqc
#SBATCH -o logs/fastqc-%j.out
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:10
#SBATCH --mem=1G
#SBATCH --array=1-8

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate RBGA_1

# Define parameters
MIN_QUAL=40
MIN_LEN=275
DOWNSAMPLING=150000
OUTPUT_DIR="05_fastQC"
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
OUTPUT_DIR="05_fastQC/$STRAIN/"$MIN_QUAL"_"$MIN_LEN"/$DOWNSAMPLING"
mkdir -p $OUTPUT_DIR

# Run fastQC
fastqc -o $OUTPUT_DIR 04_seqtk/$STRAIN/"$MIN_QUAL"_"$MIN_LEN"/$DOWNSAMPLING/outR*P_downsampled.fastq
