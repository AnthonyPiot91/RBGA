#!/bin/bash
#SBATCH -J basecalling
#SBATCH -o logs/basecalling-%j.out
#SBATCH -c 1
#SBATCH -p ibis_small
#SBATCH --time=0-05:00
#SBATCH --mem=10G
#SBATCH --gres=gpu:4

# Must run on ibis partitions because of cuda library compatibility

# Define parameters
# Change the barcode kit according to the one used for sequencing
FLOWCELL="FLO-MIN106"
echo $FLOWCELL
BARCODE_KIT="SQK-PBK004"
echo $BARCODE_KIT
FAST5_PATH=$(sed -n 5p ../config.txt)
echo $FAST5_PATH
OUTPUT_DIR="01_basecalling"

# Run guppy
/mnt/ibis/rclevesq/software/ont-guppy-v6.4.6/bin/guppy_basecaller \
	-r \
	-i $FAST5_PATH \
	-s $OUTPUT_DIR \
	--flowcell $FLOWCELL \
	--kit "$BARCODE_KIT" \
	-x cuda:0 \
	--gpu_runners_per_device 16
