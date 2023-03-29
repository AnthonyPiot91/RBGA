#!/bin/bash
#SBATCH -D ./
#SBATCH -J fastqc
#SBATCH -o logs/fastqc-%j.out
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:10
#SBATCH --mem=1G
#SBATCH --array=1-5

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate HBAW_1

# define parameters
OUTPUT_DIR="03_fastQC"

STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

# make output subdir
if [ ! -d $OUTPUT_DIR/$STRAIN ]
        then mkdir -p $OUTPUT_DIR/$STRAIN
fi

# run fastQC
fastqc -o $OUTPUT_DIR/$STRAIN 02_trimmomatic/$STRAIN/outR*P.fastq
