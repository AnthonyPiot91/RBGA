#!/bin/bash
#SBATCH -J flye
#SBATCH -o logs/flye-%j.out
#SBATCH -c 6
#SBATCH -p small
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:45
#SBATCH --mem=10G
#SBATCH --array=4

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate flye

# Define parameters

STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
THREADS=8
echo $THREADS
REPLICONS="11_replicons"
echo $REPLICONS
PARAM="2000_15"
echo $PARAM
MAPPING="asm10"
echo $MAPPING
VERIFIED="verified_minimap2"
echo $VERIFIED
INPUT_DIR="13_minimap2/$STRAIN/$VERIFIED/$REPLICONS"
echo $INPUT_DIR
FASTQ="$INPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag_unmapped.fastq"
echo $FASTQ
OUTPUT_DIR="09_flye/$STRAIN/verified/$REPLICONS/unmapped/$VERIFIED/"$PARAM"_"$MAPPING""
mkdir -p $OUTPUT_DIR

# run flye
flye 	--nano-hq \
	$FASTQ \
	--meta \
	-i 1 \
        --read-error 0.032 \
	--threads $THREADS \
	--out-dir $OUTPUT_DIR/
