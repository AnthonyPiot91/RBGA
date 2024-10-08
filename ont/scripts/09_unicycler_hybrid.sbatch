#!/bin/bash
#SBATCH -J unicycler_hybrid
#SBATCH -o logs/unicycler_hybrid-%j.out
#SBATCH -c 16
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-02:00
#SBATCH --mem=15G
#SBATCH --array=2-8

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate unicycler

# define paramaters
THREADS=16
MIN_QUAL=40
MIN_LEN=275
echo $THREADS
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
PARAM="2000_15"
echo $PARAM
LONG_READS="06_chopper/"$STRAIN"_trimmed_"$PARAM".fastq"
echo $LONG_READS
SMALL_READ_PATH="../miseq/02_trimmomatic/$STRAIN/"$MIN_QUAL"_"$MIN_LEN""
echo $SMALL_READ_PATH
SMALL_READS_R1="$SMALL_READ_PATH/outR1P.fastq"
echo $SMALL_READS_R1
SMALL_READS_R2="$SMALL_READ_PATH/outR2P.fastq"
echo $SMALL_READS_R2
OUT_DIR="09_unicycler_hybrid/$STRAIN/"$MIN_QUAL"_"$MIN_LEN"/$PARAM"
mkdir -p $OUT_DIR

# run unicycler
unicycler \
	-1 $SMALL_READS_R1 \
	-2 $SMALL_READS_R2 \
	-l $LONG_READS \
	-o $OUT_DIR \
	-t $THREADS
