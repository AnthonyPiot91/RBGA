#!/bin/bash
#SBATCH -J demultiplexing
#SBATCH -o logs/demultiplexing-%j.out
#SBATCH -c 8
#SBATCH -p ibis_small
#SBATCH --time=00-01:00
#SBATCH --mem=4G

# Must run on ibis partitions because of cuda library compatibility

# Define parameters
# Change the barcode kit according to the one used during sequencing
THREADS=16
BARCODE_KIT="SQK-PBK004"

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate HBAW_1

INPUT_DIR="01_basecalling/pass"
OUTPUT_DIR="02_demultiplexing"

# Run guppy
/mnt/ibis/rclevesq/software/ont-guppy-v6.4.6/bin/guppy_barcoder \
	-t $THREADS \
	-i $INPUT_DIR \
	-s $OUTPUT_DIR \
	--barcode_kits "$BARCODE_KIT"
