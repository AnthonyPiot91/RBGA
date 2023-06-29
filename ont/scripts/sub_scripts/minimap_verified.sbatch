#!/bin/bash
#SBATCH -J polishing
#SBATCH -o logs/polishing-%j.out
#SBATCH -c 6
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=5G
#SBATCH --array=2

# Define parameters
THREADS=12
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
FASTQ="06_chopper/"$STRAIN"_trimmed_2000_10.fastq"
echo $FASTQ
ASSEMBLY="verified"
echo $ASSEMBLY
REPLICONS="13_replicons"
echo $REPLICONS
REF="09_flye/$ASSEMBLY/$STRAIN/$REPLICONS.fasta"
echo $REF

# Load conda environment
source /mnt/ibis/rclevesq/software/miniconda3/anpio1/etc/profile.d/conda.sh

# Activate minimap2 environment -------------------------------------
conda deactivate
conda activate minimap2
echo -e "\n## Minimap2 ##\n"

# Capture file names
OUTPUT_DIR="10_minimap2/$ASSEMBLY/$REPLICONS"
mkdir -p $OUTPUT_DIR
SUFFIX="map-ont"

minimap2 --version

# Run minimap2
minimap2 -a \
         -x map-ont \
         $REF \
         $FASTQ \
         -o $OUTPUT_DIR/"$STRAIN"_"$SUFFIX".sam

# Activate samtools environment -------------------------------------
conda deactivate
conda activate samtools
echo -e "\n## SAMtools ##\n"

samtools --version | sed -n '1,3p'

# Index reference file
samtools faidx $REF

# Reflag poorly aligned reads
python scripts/sub_scripts/flag_unmapped.py \
        $REF.fai \
        $OUTPUT_DIR/"$STRAIN"_"$SUFFIX".sam \
        $OUTPUT_DIR/"$STRAIN"_"$SUFFIX"_reflag.sam
