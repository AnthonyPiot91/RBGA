#!/bin/bash
#SBATCH -J polishing
#SBATCH -o logs/polishing-%j.out
#SBATCH -c 6
#SBATCH -p ibis_small
#SBATCH --mail-type=NONE
#SBATCH --time=0-00:45
#SBATCH --mem=4G
#SBATCH --array=1

# Define parameters
THREADS=12
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
PARAM="2000_10"
echo $PARAM
FASTQ="06_chopper/"$STRAIN"_trimmed_"$PARAM".fastq"
echo $FASTQ
ASSEMBLY="unicycler_hybrid_all_clean_reads"
echo $ASSEMBLY
REF="09_unicycler_hybrid/$STRAIN/all_clean_reads/assembly.fasta"
echo $REF
MAPPING="asm10"
echo $MAPPING
RACON_ERR_THR="0.1" #0.03 or 0.1
echo $RACON_ERR_THR

# Load conda environment
source /mnt/ibis/rclevesq/software/miniconda3/anpio1/etc/profile.d/conda.sh

# Activate bbmap environment ----------------------------------------
conda activate bbmap
echo -e "\n## Tadpole ##\n"

OUTPUT_DIR="10_tadpole/$STRAIN/$ASSEMBLY"
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
                out=$OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_tadpole1.fa \
                mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=$OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_tadpole1.fa \
                extra=$FASTQ \
                out=$OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_tadpole2.fa \
                mode=extend

# Activate minimap2 environment -------------------------------------
conda deactivate
conda activate minimap2
echo -e "\n## Minimap2 ##\n"

# Capture file names
REF="$OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_tadpole2.fa"
echo $REF
OUTPUT_DIR="11_minimap2/$STRAIN/$ASSEMBLY"
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

# Activate racon environment ----------------------------------------
conda deactivate
conda activate racon
echo -e "\n## Racon ##\n"

# Capture file names
SAM="$OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_recall_reflag.sam"
echo $SAM
OUTPUT_DIR="12_racon/$STRAIN/$ASSEMBLY"
mkdir -p $OUTPUT_DIR

racon --version

# Run racon
racon   $FASTQ \
        $SAM \
        $REF \
        -t $THREADS \
        --error-threshold $RACON_ERR_THR \
        --no-trimming \
        > $OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_racon.fa

# Activate minimap2 environment -------------------------------------
conda deactivate
conda activate minimap2
echo -e "\n## Minimap2 ##\n"

# Capture file names
REF="$OUTPUT_DIR/"$STRAIN"_"$PARAM"_"$MAPPING"_racon.fa"
echo $REF
OUTPUT_DIR="13_minimap2/$STRAIN/$ASSEMBLY"
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
