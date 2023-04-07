#!/bin/bash
#SBATCH -J canu
#SBATCH -o logs/canu-%j.out
#SBATCH -c 6
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-05:00
#SBATCH --mem=10G
#SBATCH --array=1-5

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate canu

# Define parameters
THREADS=16
GENOME_SIZE="1.5m"
INPUT_DIR="06_chopper"
OUTPUT_DIR="09_canu/with_chrom"

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
BARCODE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f3)
echo $BARCODE
INPUT_FILE_PATH="$INPUT_DIR/"$STRAIN"*.fastq"
echo $INPUT_FILE_PATH

# Get canu version
canu --version

# Run canu
canu 	-p $STRAIN \
	-d $OUTPUT_DIR/"$STRAIN" \
	genomeSize=$GENOME_SIZE \
	maxThreads=$THREADS \
	useGrid=false \
	-nanopore $INPUT_FILE_PATH
