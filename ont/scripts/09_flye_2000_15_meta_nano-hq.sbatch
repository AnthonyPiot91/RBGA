#!/bin/bash
#SBATCH -J flye
#SBATCH -o logs/flye-%j.out
#SBATCH -c 6
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-01:30
#SBATCH --mem=10G
#SBATCH --array=1-8

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate flye

STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

# define paramaters
THREADS=8
INPUT_DIR="06_chopper"
OUTPUT_DIR="09_flye/$STRAIN/2000_15_meta_nano-hq"
GENOME_SIZE=1500000
echo $GENOME_SIZE
INPUT_FILE_PATH="$INPUT_DIR/"$STRAIN"_trimmed_2000_15.fastq"
echo $INPUT_FILE_PATH

if [ ! -d $OUTPUT_DIR ]
    then mkdir -p $OUTPUT_DIR
fi

# run flye
flye 	--nano-hq \
	$INPUT_FILE_PATH \
	-i 2 \
	--meta \
	--genome-size $GENOME_SIZE \
        --read-error 0.03 \
	--threads $THREADS \
	--out-dir $OUTPUT_DIR/
