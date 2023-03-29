#!/bin/bash
#SBATCH -D ./
#SBATCH -J A5-miseq
#SBATCH -o logs/A5-miseq-%j.out
#SBATCH -c 8
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-05:00
#SBATCH --mem=5G
#SBATCH --array=1-5

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate HBAW_1

# define parameters
THREADS=8
OUTPUT_DIR="04_A5-miseq"

STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

MISEQ_ID=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f2)
echo $MISEQ_ID

RAW_SMALL_READ_PATH=$(sed -n 2p ../config.txt | sed 's|/$||')
echo $RAW_SMALL_READ_PATH

# make output subdir
if [ ! -d $OUTPUT_DIR/$STRAIN ]
	then mkdir -p $OUTPUT_DIR/$STRAIN
fi

cd $OUTPUT_DIR/$STRAIN

# run A5-miseq
a5_pipeline.pl 	--threads $THREADS \
		$RAW_SMALL_READ_PATH/*"$MISEQ_ID"*R1*.fastq.gz \
		$RAW_SMALL_READ_PATH/*"$MISEQ_ID"*R2*.fastq.gz \
		$STRAIN
