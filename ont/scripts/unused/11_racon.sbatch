#!/bin/bash
#SBATCH -J racon
#SBATCH -o logs/racon-%j.out
#SBATCH -c 8
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=5G
#SBATCH --array=2

# Load conda environment
source /mnt/ibis/rclevesq/software/miniconda3/anpio1/etc/profile.d/conda.sh
conda activate racon

# Define paramaters
THREADS=16

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
REF="09_flye/2000_15_meta_nano-hq/$STRAIN/assembly.fasta"
echo $REF
SAM="10_minimap2/flye_2000_15_meta_nano-hq/$STRAIN/"$STRAIN"_flye_2000_10_meta_nano-hq_map-ont_reflag.sam"
echo $SAM

racon --version

# Run racon
racon	06_chopper/"$STRAIN"_trimmed_2000_10.fastq \
	$SAM \
	$REF \
	-t $THREADS \
	> 11_racon/"$STRAIN"_racon.fa
