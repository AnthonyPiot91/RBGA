#!/bin/bash
#SBATCH -J demultiplexing
#SBATCH -o logs/demultiplexing-%j.out
#SBATCH -c 8
#SBATCH -p ibis_small
#SBATCH --time=00-01:00
#SBATCH --mem=4G

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate RBGA_1

# Define parameters
# Change the barcode kit according to the one used during sequencing
THREADS=16
BARCODE_KIT="EXP-NBD104 EXP-NBD114"
INPUT_DIR="01_basecalling/pass"
OUTPUT_DIR="02_demultiplexing"

# Run guppy
/mnt/ibis/rclevesq/software/ont-guppy-v6.4.6/bin/guppy_barcoder \
	-t $THREADS \
	-i $INPUT_DIR \
	-s $OUTPUT_DIR \
	--barcode_kits "$BARCODE_KIT" \
	--detect_adapter \
        --detect_barcodes \
        --detect_mid_strand_adapter \
        --detect_mid_strand_barcodes \
        --detect_primer \
        --enable_trim_barcodes \
        --trim_adapters \
        --trim_primers
