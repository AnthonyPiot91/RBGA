#!/bin/bash
#SBATCH -J tadpole
#SBATCH -o logs/tadpole-%j.out
#SBATCH -c 6
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:30
#SBATCH --mem=6G
#SBATCH --array=2

# Load conda environment
source /mnt/ibis/rclevesq/software/miniconda3/anpio1/etc/profile.d/conda.sh
conda activate bbmap

# Define paramaters
THREADS=12

# Capture fastq file name
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
REF="09_flye/2000_15_meta_nano-hq/$STRAIN/assembly.fasta"
echo $REF

tadpole.sh --version

# Run 5 round of contig extension with tadpole
tadpole.sh	threads=$THREADS \
		extendleft=1000 \
		extendright=1000 \
		trimcircular=t \
		overwrite=true \
		in=$REF \
		extra=06_chopper/"$STRAIN"_trimmed_2000_10.fastq \
		out=tadpole/"$STRAIN"_tadpole1.fa \
		mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=tadpole/"$STRAIN"_tadpole1.fa \
                extra=06_chopper/"$STRAIN"_trimmed_2000_10.fastq \
                out=tadpole/"$STRAIN"_tadpole2.fa \
                mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=tadpole/"$STRAIN"_tadpole2.fa \
                extra=06_chopper/"$STRAIN"_trimmed_2000_10.fastq \
                out=tadpole/"$STRAIN"_tadpole3.fa \
                mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=tadpole/"$STRAIN"_tadpole3.fa \
                extra=06_chopper/"$STRAIN"_trimmed_2000_10.fastq \
                out=tadpole/"$STRAIN"_tadpole4.fa \
                mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=tadpole/"$STRAIN"_tadpole4.fa \
                extra=06_chopper/"$STRAIN"_trimmed_2000_10.fastq \
                out=tadpole/"$STRAIN"_tadpole5.fa \
                mode=extend
