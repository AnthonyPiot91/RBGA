#!/bin/bash
#SBATCH -J mafft
#SBATCH -o logs/mafft-%j.out
#SBATCH -c 2
#SBATCH -p ibis_small
#SBATCH --time=0-01:00
#SBATCH --mem=1G

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate mafft

mafft 	--thread 8 \
	--reorder \
	mafft/cp32_long_reads.fasta \
	> mafft/cp32_long_reads_aligned.fasta

