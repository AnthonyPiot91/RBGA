#!/bin/bash
#SBATCH -J cat
#SBATCH -o logs/cat-%j.out
#SBATCH -c 6
#SBATCH -p batch
#SBATCH --mail-type=NONE
#SBATCH --time=00-01:00
#SBATCH --mem=1G

# Define parameters
INPUT_DIR="02_demultiplexing"
echo $INPUT_DIR
OUTPUT_DIR="03_concatenate"
echo $OUTPUT_DIR

# Concatenate barcoding files
ls $INPUT_DIR/ | grep "barcode[0-9][0-9]" | sort -Vu | \
while read barcode
        do echo "cat $INPUT_DIR/$barcode/*.fastq > $OUTPUT_DIR/$barcode.merged.fastq"
done | parallel -k -j6
