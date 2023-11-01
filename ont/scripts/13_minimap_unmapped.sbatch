#!/bin/bash
#SBATCH -J minimap2
#SBATCH -o logs/minimap2-%j.out
#SBATCH -c 10
#SBATCH -p small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:15
#SBATCH --mem=4G
#SBATCH --array=1

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh

# Define parameters
THREADS=10
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
PARAM="2000_13"
echo $PARAM
MAPPING="asm5"
echo $MAPPING
ASSEMBLY="verified"
echo $ASSEMBLY
ORIGINAL_POLISHING=""$ASSEMBLY"_minimap2"
echo $ORIGINAL_POLISHING
POLISHING=""$ASSEMBLY"_minimap2"
echo $POLISHING
REPLICONS="12_replicons"
echo $REPLICONS
REF="09_flye/$STRAIN/$ASSEMBLY/$REPLICONS/unmapped/$ORIGINAL_POLISHING/"$PARAM"_"$MAPPING"/assembly.fasta"
echo $REF
PREFIX=""$STRAIN"_"$PARAM"_"$MAPPING""
echo $PREFIX
FASTQ="13_minimap2/$STRAIN/$ORIGINAL_POLISHING/$REPLICONS/"$PREFIX"_recall_reflag_unmapped.fastq"
echo $FASTQ
OUTPUT_DIR="13_minimap2/$STRAIN/$POLISHING/$REPLICONS/unmapped"
mkdir -p $OUTPUT_DIR

echo -e "\n## Minimap2 ##\n" #--------------------------------------------------

# Activate minimap2 environment
conda activate minimap2

# Run minimap2 script
bash scripts/sub_scripts/minimap.sh \
    $MAPPING \
    $REF \
    $FASTQ \
    $OUTPUT_DIR/$PREFIX

echo -e "\n## SAMtools ##\n" #--------------------------------------------------

# Activate samtools environment
conda activate samtools

# Run samtools script
bash scripts/sub_scripts/samtools.sh \
    $MAPPING \
    $REF \
    $FASTQ \
    $OUTPUT_DIR/$PREFIX \
    $THREADS
