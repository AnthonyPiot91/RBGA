#!/bin/bash
#SBATCH -J flye
#SBATCH -o logs/flye-%j.out
#SBATCH -c 6
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=00-00:45
#SBATCH --mem=10G
#SBATCH --array=5

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate flye

STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

# define paramaters
THREADS=8
echo $THREADS
REPLICONS="18_replicons"
echo $REPLICONS
PARAM="2000_10"
echo $PARAM
MAPPING="asm10"
echo $MAPPING
FASTQ="13_minimap2/$STRAIN/verified_minimap2_only/$REPLICONS/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag_unmapped.fastq"
echo $FASTQ
OUTPUT_DIR="09_flye/$STRAIN/verified/$REPLICONS/unmapped/"$PARAM"_"$MAPPING""
mkdir -p $OUTPUT_DIR

# run flye
flye 	--nano-hq \
	$FASTQ \
	--meta \
	-i 1 \
        --read-error 0.1 \
	--threads $THREADS \
	--out-dir $OUTPUT_DIR/
