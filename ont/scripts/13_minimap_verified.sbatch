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
REPLICONS="12_replicons"
echo $REPLICONS
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
PARAM="2000_15"
echo $PARAM
MAPPING="asm10"
echo $MAPPING
FASTQ="06_chopper/"$STRAIN"_trimmed_"$PARAM".fastq"
echo $FASTQ
ASSEMBLY="verified"
echo $ASSEMBLY
POLISHING=""$ASSEMBLY"_minimap2"
echo $POLISHING
REF="10_verified/$STRAIN/$REPLICONS/$REPLICONS.fasta"
echo $REF
PREFIX=""$STRAIN"_"$PARAM"_"$MAPPING""
echo $PREFIX
OUTPUT_DIR="13_minimap2/$STRAIN/$POLISHING/$REPLICONS"
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
