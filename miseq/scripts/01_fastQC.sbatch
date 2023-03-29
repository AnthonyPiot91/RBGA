#!/bin/bash
#SBATCH -D ./
#SBATCH -J fastqc
#SBATCH -o logs/fastqc-%j.out
#SBATCH -c 4
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:10
#SBATCH --mem=5G

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate HBAW_1

# define parameters
THREADS=4
OUTPUT_DIR="01_fastQC"

# capture raw Miseq fastq file names
RAW_SMALL_READ_PATH=$(sed -n 2p ../config.txt | sed 's|/$||')
echo $RAW_SMALL_READ_PATH

# run fastQC
fastqc -t $THREADS -o $OUTPUT_DIR/ $RAW_SMALL_READ_PATH/*.fastq.gz
