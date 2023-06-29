#!/bin/bash
#SBATCH -J flye
#SBATCH -o logs/flye-%j.out
#SBATCH -c 6
#SBATCH -p batch
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
REPLICONS="11_replicons"
FASTQ="10_minimap2/$STRAIN/verified/$REPLICONS/"$STRAIN"_map-ont_reflag_unmapped.fastq"
echo $FASTQ
OUTPUT_DIR="09_flye/$STRAIN/verified/$REPLICONS/unmapped"

# run flye
flye 	--nano-hq \
	$FASTQ \
	-i 2 \
	--meta \
        --read-error 0.06 \
	--threads $THREADS \
	--out-dir $OUTPUT_DIR/
