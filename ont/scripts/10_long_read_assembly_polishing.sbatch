#!/bin/bash
#SBATCH -J polishing
#SBATCH -o logs/polishing-%j.out
#SBATCH -c 6
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:45
#SBATCH --mem=4G
#SBATCH --array=5

# Define parameters
THREADS=12
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
FASTQ="06_chopper/"$STRAIN"_trimmed_2000_10.fastq"
echo $FASTQ
ASSEMBLY="2000_10_meta_nano-hq"
echo $ASSEMBLY
REPLICONS="assembly"
echo $REPLICONS
REF="09_flye/$STRAIN/$ASSEMBLY/$REPLICONS.fasta"
echo $REF

# Load conda environment
source /mnt/ibis/rclevesq/software/miniconda3/anpio1/etc/profile.d/conda.sh

# Activate bbmap environment ----------------------------------------
conda activate bbmap
echo -e "\n## Tadpole ##\n"

OUTPUT_DIR="10_tadpole/$STRAIN/$ASSEMBLY/$REPLICONS"
mkdir -p $OUTPUT_DIR

tadpole.sh --version

# Run 5 round of contig extension with tadpole
tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=$REF \
                extra=$FASTQ \
                out=$OUTPUT_DIR/"$STRAIN"_tadpole1.fa \
                mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=$OUTPUT_DIR/"$STRAIN"_tadpole1.fa \
                extra=$FASTQ \
                out=$OUTPUT_DIR/"$STRAIN"_tadpole2.fa \
                mode=extend

# Activate minimap2 environment -------------------------------------
conda deactivate
conda activate minimap2
echo -e "\n## Minimap2 ##\n"

# Capture file names
REF="$OUTPUT_DIR/"$STRAIN"_tadpole2.fa"
echo $REF
OUTPUT_DIR="11_minimap2/$STRAIN/$ASSEMBLY/$REPLICONS"
mkdir -p $OUTPUT_DIR
SUFFIX="map-ont"

minimap2 --version

# Run minimap2
minimap2 -a \
         -x map-ont \
         $REF \
         $FASTQ \
         -o $OUTPUT_DIR/"$STRAIN"_"$SUFFIX".sam

# Activate samtools environment -------------------------------------
conda deactivate
conda activate samtools
echo -e "\n## SAMtools ##\n"

samtools --version | sed -n '1,3p'

# Index reference file
samtools faidx $REF

# Reflag poorly aligned reads
python scripts/sub_scripts/flag_unmapped.py \
        $REF.fai \
        $OUTPUT_DIR/"$STRAIN"_"$SUFFIX".sam \
        $OUTPUT_DIR/"$STRAIN"_"$SUFFIX"_reflag.sam

# Activate racon environment ----------------------------------------
conda deactivate
conda activate racon
echo -e "\n## Racon ##\n"

# Capture file names
SAM="$OUTPUT_DIR/"$STRAIN"_"$SUFFIX"_reflag.sam"
echo $SAM
OUTPUT_DIR="12_racon/$STRAIN/$ASSEMBLY/$REPLICONS"
mkdir -p $OUTPUT_DIR

racon --version

# Run racon
racon   $FASTQ \
        $SAM \
        $REF \
        -t $THREADS \
        --error-threshold 0.1 \
        --no-trimming \
        > $OUTPUT_DIR/"$STRAIN"_racon.fa

# Activate minimap2 environment -------------------------------------
conda deactivate
conda activate minimap2
echo -e "\n## Minimap2 ##\n"

# Capture file names
REF="$OUTPUT_DIR/"$STRAIN"_racon.fa"
echo $REF
OUTPUT_DIR="13_minimap2/$STRAIN/$ASSEMBLY/$REPLICONS"
mkdir -p $OUTPUT_DIR
SUFFIX="map-ont"

minimap2 --version

# Run minimap2
minimap2 -a \
	 -x map-ont \
	 $REF \
	 $FASTQ \
	 -o $OUTPUT_DIR/"$STRAIN"_"$SUFFIX".sam

# Activate samtools environment -------------------------------------
conda deactivate
conda activate samtools
echo -e "\n## SAMtools ##\n"

samtools --version | sed -n '1,3p'

# Index reference file
samtools faidx $REF

# Reflag poorly aligned reads
python scripts/sub_scripts/flag_unmapped.py \
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
