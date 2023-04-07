#!/bin/bash
#SBATCH -J stats
#SBATCH -o logs/stats-%j.out
#SBATCH -c 1
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=500M
#SBATCH --array=1-5

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh

# define paramaters
INPUT_DIR="10_minimap2"
OUTPUT_DIR="14_statistics"

# capture info on input files
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

echo -e "Method\tPercentage of unmapped reads" > $OUTPUT_DIR/$STRAIN"_unmapped_reads.txt"

for dir in "with_chrom" "polypolish_canu" "quickmerge_canu" "quickmerge_polypolish_canu"
	do
	RAW_READS=$(samtools view -c -f 0 $INPUT_DIR/$dir/$STRAIN*.sam)
	UNMAPPED_READS=$(samtools view -c -f 4 $INPUT_DIR/$dir/$STRAIN*.sam)
	percentage=$(echo "scale=3; ($UNMAPPED_READS * 100) / $RAW_READS" | \
	bc | sed 's/$/%/ ; s/^\./0./')
	echo -e $dir"\t"$percentage

done >> $OUTPUT_DIR/$STRAIN"_unmapped_reads.txt"
