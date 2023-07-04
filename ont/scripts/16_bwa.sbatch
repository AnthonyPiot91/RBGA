#!/bin/bash
#SBATCH -J bwa
#SBATCH -o logs/bwa-%j.out
#SBATCH -c 8
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=10G
#SBATCH --array=4

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate RBGA_1

# Define paramaters
THREADS=10
OUTPUT_DIR="16_bwa"
SMALL_READ_PATH="../miseq/02_trimmomatic"

# Capture strain an barcode name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
BARCODE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f3)
echo $BARCODE

# Capture long read fasta assembly file name
REF="15_polypolish/$STRAIN/"$STRAIN"_polished.fasta"
echo $REF

# Capture small read fastq name
SMALL_READS_R1=$(ls -1 $SMALL_READ_PATH/$STRAIN/outR1P.fastq)
SMALL_READS_R2=$(ls -1 $SMALL_READ_PATH/$STRAIN/outR2P.fastq)
echo $SMALL_READS_R1
echo $SMALL_READS_R2

# Run bwa-mem
bwa index $REF
bwa mem -t $THREADS -a $REF $SMALL_READS_R1 $SMALL_READS_R2> $OUTPUT_DIR/"$STRAIN"_polished.sam

# Activate samtools
conda deactivate
conda activate samtools

# Sort sam file
samtools sort 	-@ $THREADS \
		-O SAM \
		$OUTPUT_DIR/"$STRAIN"_polished.sam \
		-o $OUTPUT_DIR/"$STRAIN"_polished_sorted.sam

rm $OUTPUT_DIR/"$STRAIN"_polished.sam
