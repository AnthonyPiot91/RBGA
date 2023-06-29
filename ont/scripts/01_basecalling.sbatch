#!/bin/bash
#SBATCH -J basecalling
#SBATCH -o logs/basecalling-%j.out
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --time=0-10:00
#SBATCH --mem=5G
#SBATCH --gres=gpu:4

FAST5_PATH=$(sed -n 5p ../config.txt)
echo $FAST5_PATH

OUTPUT_DIR="01_basecalling_sup"
FLOWCELL="FLO-MIN111"
echo $FLOWCELL
BARCODE_KIT="SQK-LSK109"
echo $BARCODE_KIT

# run guppy
/mnt/ibis/rclevesq/software/ont-guppy-v6.4.6/bin/guppy_basecaller \
	-r \
	-i $FAST5_PATH \
	-s $OUTPUT_DIR \
	-c  dna_r10.3_450bps_sup.cfg \
	-x auto \
	--gpu_runners_per_device 16
