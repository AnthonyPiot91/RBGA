#!/bin/bash
#SBATCH -J stats
#SBATCH -o logs/stats-%j.out
#SBATCH -c 1
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:15
#SBATCH --mem=500M
#SBATCH --array=2,3

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh

# define paramaters
INPUT_DIR="15_minimap2"
OUTPUT_DIR="19_statistics"

# capture info on input files
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

REPLICONS=$(cat $OUTPUT_DIR/$STRAIN)

echo -e "Strain\tRaw reads\tUnmapped reads\tPercentage of unmapped reads" > $OUTPUT_DIR/$STRAIN"_unmapped_reads.txt"

RAW_READS=$(samtools view -c -f 0 $INPUT_DIR/$STRAIN/verified/$REPLICONS/"$STRAIN"*reflag.sam)
UNMAPPED_READS=$(samtools view -c -f 4 $INPUT_DIR/$STRAIN/verified/$REPLICONS/"$STRAIN"*reflag.sam)
percentage=$(echo "scale=3; ($UNMAPPED_READS * 100) / $RAW_READS" | bc | sed 's/$/%/ ; s/^\./0./')

echo -e "$STRAIN\t$RAW_READS\t$UNMAPPED_READS\t$percentage" >> $OUTPUT_DIR/$STRAIN"_unmapped_reads.txt"
