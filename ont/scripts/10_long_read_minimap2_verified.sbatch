#!/bin/bash
#SBATCH -J polishing_verified
#SBATCH -o logs/polishing_verified-%j.out
#SBATCH -c 6
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:15
#SBATCH --mem=4G
#SBATCH --array=2

# Define parameters
THREADS=12
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
PARAM="2000_10"
echo $PARAM
FASTQ="06_chopper/"$STRAIN"_trimmed_"$PARAM".fastq"
echo $FASTQ
ASSEMBLY="verified"
echo $ASSEMBLY
REPLICONS="19_replicons"
echo $REPLICONS
REF="09_flye/$STRAIN/$ASSEMBLY/$REPLICONS/$REPLICONS.fasta"
echo $REF
MAPPING="asm10"
echo $MAPPING
RACON_ERR_THR="0.1" #0.03 or 0.1
echo $RACON_ERR_THR

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh

# Activate minimap2 environment -------------------------------------
conda deactivate
conda activate minimap2
echo -e "\n## Minimap2 ##\n"

# Capture file names
OUTPUT_DIR="13_minimap2/$STRAIN/"$ASSEMBLY"_minimap2_only/$REPLICONS"
mkdir -p $OUTPUT_DIR

minimap2 --version

# Run minimap2
minimap2 -a \
	 -x $MAPPING \
	 $REF \
	 $FASTQ \
	 -o $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING".sam

# Activate samtools environment -------------------------------------
conda deactivate
conda activate samtools
echo -e "\n## SAMtools ##\n"

samtools --version | sed -n '1,3p'

# Index reference file
samtools faidx $REF

# Convert SAM to BAM
samtools view -b -@10 \
                 $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING".sam \
                 > $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING".bam

# Recall mismatch
samtools calmd -@10 --output-fmt SAM \
                    $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING".bam \
                    $REF \
                    > $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall.sam

# Reflag poorly aligned reads
python scripts/sub_scripts/flag_unmapped.py \
        $REF.fai \
        $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall.sam \
        $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag.sam

rm $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING".bam
rm $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING".sam

# Extract unmapped reads with samtools
samtools view -bf 4 \
                $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag.sam \
                > $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag_unmapped.bam

# Convert bam file to fastq
samtools fastq $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag_unmapped.bam \
             > $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag_unmapped.fastq

rm $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag_unmapped.bam
