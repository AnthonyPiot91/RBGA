#!/bin/bash
#SBATCH -D ./
#SBATCH -J trimmomatic
#SBATCH -o logs/trimmomatic-%j.out
#SBATCH -c 4
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:30
#SBATCH --mem=1G
#SBATCH --array=1-8

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate RBGA_1

# define paramaters
THREADS=8
MIN_QUAL=40
MIN_LEN=275
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
MISEQ_ID=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f2)
echo $MISEQ_ID
RAW_SMALL_READ_PATH=$(sed -n 2p ../config.txt | sed 's|/$||')
echo $RAW_SMALL_READ_PATH
OUTPUT_DIR="02_trimmomatic/$STRAIN/"$MIN_QUAL"_"$MIN_LEN""
mkdir -p $OUTPUT_DIR

# run Trimmomatic
trimmomatic PE -threads $THREADS \
			$RAW_SMALL_READ_PATH/$MISEQ_ID/"$MISEQ_ID"*R1*.fastq.gz \
			$RAW_SMALL_READ_PATH/$MISEQ_ID/"$MISEQ_ID"*R2*.fastq.gz \
			$OUTPUT_DIR/outR1P.fastq \
			$OUTPUT_DIR/outR1U.fastq \
			$OUTPUT_DIR/outR2P.fastq \
			$OUTPUT_DIR/outR2U.fastq \
			ILLUMINACLIP:TruSeq3-PE.fa:4:20:10 \
			LEADING:$MIN_QUAL TRAILING:$MIN_QUAL \
			SLIDINGWINDOW:8:$MIN_QUAL MINLEN:$MIN_LEN
