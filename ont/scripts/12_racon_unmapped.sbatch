#!/bin/bash
#SBATCH -J racon
#SBATCH -o logs/racon-%j.out
#SBATCH -c 10
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:15
#SBATCH --mem=4G
#SBATCH --array=2

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh

# Define parameters
THREADS=10
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
PARAM="2000_15"
echo $PARAM
MAPPING="asm5"
echo $MAPPING
ASSEMBLY="verified"
echo $ASSEMBLY
ORIGINAL_POLISHING=""$ASSEMBLY"_racon"
echo $ORIGINAL_POLISHING
POLISHING=""$ASSEMBLY"_racon"
echo $POLISHING
REPLICONS="15_replicons"
echo $REPLICONS
REF="09_flye/$STRAIN/$ASSEMBLY/$REPLICONS/unmapped/$ORIGINAL_POLISHING/"$PARAM"_"$MAPPING"/assembly.fasta"
echo $REF
PREFIX=""$STRAIN"_"$PARAM"_"$MAPPING""
echo $PREFIX
FASTQ="13_minimap2/$STRAIN/$ORIGINAL_POLISHING/$REPLICONS/"$PREFIX"_recall_reflag_unmapped.fastq"
echo $FASTQ
OUTPUT_DIR="11_minimap2/$STRAIN/$POLISHING/$REPLICONS/unmapped"
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

echo -e "\n## Racon ##\n" #-----------------------------------------------------

# Activate racon environment
conda activate racon

# Define parameters
RACON_ERR_THR="0.032"
echo $RACON_ERR_THR
SAM="$OUTPUT_DIR/"$PREFIX"_recall_reflag.sam"
echo $SAM
OUTPUT_DIR="12_racon/$STRAIN/$POLISHING/$REPLICONS/unmapped"
mkdir -p $OUTPUT_DIR

# Run racon script
bash scripts/sub_scripts/racon.sh \
    $MAPPING \
    $REF \
    $FASTQ \
    $OUTPUT_DIR/$PREFIX \
    $THREADS \
    $RACON_ERR_THR \
    $SAM

echo -e "\n## Minimap2 ##\n" #--------------------------------------------------

# Activate minimap2 environment
conda activate minimap2

# Define parameters
REF="$OUTPUT_DIR/"$PREFIX"_racon.fa"
echo $REF
OUTPUT_DIR="13_minimap2/$STRAIN/$POLISHING/$REPLICONS/unmapped"
mkdir -p $OUTPUT_DIR

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
