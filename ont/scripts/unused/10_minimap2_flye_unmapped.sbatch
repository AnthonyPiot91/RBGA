#!/bin/bash
#SBATCH -J minimap2
#SBATCH -o logs/minimap2-%j.out
#SBATCH -c 6
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=3G
#SBATCH --array=1

# Load conda environment
source /mnt/ibis/rclevesq/software/miniconda3/anpio1/etc/profile.d/conda.sh
conda activate minimap2

# Define paramaters
THREADS=12
INPUT_DIR_FASTQ="10_minimap2/flye_2000_15_meta_nano-hq"
OUTPUT_DIR="10_minimap2/flye_2000_15_meta_nano-hq"

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
REF="09_flye/2000_15_meta_nano-hq/$STRAIN/unmapped/assembly.fasta"
echo $REF
SUFFIX="flye_2000_10_18_verified_unmapped_map-ont"

minimap2 --version

# Run minimap2
minimap2 -a \
	 -x map-ont \
	 $REF \
	 $INPUT_DIR_FASTQ/"$STRAIN"_flye_2000_10_18_verified_map-ont_reflag_unmapped.fastq \
	 -o $OUTPUT_DIR/"$STRAIN"_"$SUFFIX".sam

# Activate samtools environment
conda deactivate
conda activate samtools

samtools --version | sed -n '1,3p'

# Index reference file
samtools faidx $REF

# Reflag poorly aligned reads
python scripts/flag_unmapped.py \
	$REF.fai \
        $OUTPUT_DIR/"$STRAIN"_"$SUFFIX".sam \
        $OUTPUT_DIR/"$STRAIN"_"$SUFFIX"_reflag.sam

# Extract unmapped reads with samtools
samtools view -bf 4 \
                $OUTPUT_DIR/"$STRAIN"_"$SUFFIX"_reflag.sam \
                > $OUTPUT_DIR/"$STRAIN"_"$SUFFIX"_reflag_unmapped.bam

# Convert bam file to fastq
samtools fastq $OUTPUT_DIR/"$STRAIN"_"$SUFFIX"_reflag_unmapped.bam \
             > $OUTPUT_DIR/"$STRAIN"_"$SUFFIX"_reflag_unmapped.fastq
