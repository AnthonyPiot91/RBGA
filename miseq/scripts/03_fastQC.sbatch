#!/bin/bash
#SBATCH -D ./
#SBATCH -J fastqc
#SBATCH -o logs/fastqc-%j.out
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:10
#SBATCH --mem=1G
#SBATCH --array=1-8

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate RBGA_1

# define parameters
MIN_QUAL=40
MIN_LEN=275
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
INPUT_DIR="02_trimmomatic/$STRAIN/"$MIN_QUAL"_"$MIN_LEN""
echo $INPUT_DIR
OUTPUT_DIR="03_fastQC/$STRAIN/"$MIN_QUAL"_"$MIN_LEN""
mkdir -p $OUTPUT_DIR

# run fastQC
fastqc -o $OUTPUT_DIR $INPUT_DIR/outR*P.fastq
