#!/bin/bash
#SBATCH -J polypolish
#SBATCH -o logs/polypolish-%j.out
#SBATCH -c 5
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-01:00
#SBATCH --mem=5G
#SBATCH --array=1-6

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate HBAW_1

# Define paramaters
INPUT_DIR1="09_canu"
INPUT_DIR2="11_bwa/canu"
OUTPUT_DIR="12_polypolish/canu"

# Capture info on input files
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
BARCODE=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f3)
echo $BARCODE
LONG_READS="$INPUT_DIR1/$STRAIN/assembly.fasta"
echo $LONG_READS

# Filter the alignments
polypolish_insert_filter.py \
	--in1 $INPUT_DIR2/"$STRAIN"_alignments_1.sam \
	--in2 $INPUT_DIR2/"$STRAIN"_alignments_2.sam \
	--out1 $OUTPUT_DIR/"$STRAIN"_filtered_1.sam \
	--out2 $OUTPUT_DIR/"$STRAIN"_filtered_2.sam

# Run polypolish
polypolish --debug $OUTPUT_DIR/"$STRAIN"_polished.tsv \
	   $LONG_READS \
	   $OUTPUT_DIR/"$STRAIN"_filtered_1.sam \
	   $OUTPUT_DIR/"$STRAIN"_filtered_2.sam \
	   > $OUTPUT_DIR/"$STRAIN"_polished.fasta
