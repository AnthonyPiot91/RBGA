#!/bin/bash
#SBATCH -J quickmerge
#SBATCH -o logs/quickmerge-%j.out
#SBATCH -c 1
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=1G
#SBATCH --array=1-5

# load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate quickmerge

# define paramaters
INPUT_DIR1="09_canu/with_chrom"
INPUT_DIR2="09_canu/without_chrom"
OUTPUT_DIR="13_quickmerge/canu"

# capture info on input files
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN

cd $OUTPUT_DIR

# Run quickmerge wrapper
merge_wrapper.py -pre $STRAIN \
		 ../../$INPUT_DIR1/$STRAIN/$STRAIN.contigs.fasta \
		 ../../$INPUT_DIR2/$STRAIN/$STRAIN.contigs.fasta
