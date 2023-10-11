#!/bin/bash
#SBATCH -J seqtk
#SBATCH -o logs/seqtk-%j.out
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:10
#SBATCH --mem=1G
#SBATCH --array=1-8

# Define parameters
MIN_QUAL=40
MIN_LEN=275
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
SMALL_READS_R1="02_trimmomatic/$STRAIN/"$MIN_QUAL"_"$MIN_LEN"/outR1P.fastq"
echo $SMALL_READS_R1
SMALL_READS_R2="02_trimmomatic/$STRAIN/"$MIN_QUAL"_"$MIN_LEN"/outR2P.fastq"
echo $SMALL_READS_R2
DOWNSAMPLING=75000
echo $DOWNSAMPLING
OUTPUT_DIR="04_seqtk/$STRAIN/"$MIN_QUAL"_"$MIN_LEN"/$DOWNSAMPLING"
mkdir -p $OUTPUT_DIR

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh

# Activate seqtk environment -------------------------------------
conda activate seqtk

seqtk sample -s100 $SMALL_READS_R1 $DOWNSAMPLING > $OUTPUT_DIR/outR1P_downsampled.fastq
seqtk sample -s100 $SMALL_READS_R2 $DOWNSAMPLING > $OUTPUT_DIR/outR2P_downsampled.fastq
