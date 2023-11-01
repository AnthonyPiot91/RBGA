#!/bin/bash
#SBATCH -J flye
#SBATCH -o logs/flye-%j.out
#SBATCH -c 6
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-01:30
#SBATCH --mem=10G
#SBATCH --array=2-8

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate flye

STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

# define paramaters
THREADS=8
INPUT_DIR="06_chopper"
GENOME_SIZE=1500000
PARAM="4000_15"
echo $GENOME_SIZE
INPUT_FILE_PATH="$INPUT_DIR/"$STRAIN"_trimmed_"$PARAM".fastq"
echo $INPUT_FILE_PATH
OUTPUT_DIR="09_flye/$STRAIN/"$PARAM"_meta_nano-hq"
mkdir -p $OUTPUT_DIR

# run flye
flye 	--nano-hq \
	$INPUT_FILE_PATH \
	-i 2 \
	--meta \
	--genome-size $GENOME_SIZE \
        --read-error 0.032 \
	--threads $THREADS \
	--out-dir $OUTPUT_DIR/
