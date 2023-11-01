# Define parameters
MAPPING=$1
REF=$2
FASTQ=$3
OUTPUT=$4
THREADS=$5

# Print samtools version
samtools --version | sed -n '1,3p'

# Index reference file
samtools faidx $REF

# Convert SAM to BAM
samtools view -b \
              -@ $THREADS \
              $OUTPUT.sam \
              > $OUTPUT.bam

# Recall mismatch
samtools calmd -@ $THREADS \
                --output-fmt SAM \
                $OUTPUT.bam \
                $REF \
                > "$OUTPUT"_recall.sam

# Reflag poorly aligned reads
python scripts/sub_scripts/flag_unmapped.py \
        $REF.fai \
        "$OUTPUT"_recall.sam \
        "$OUTPUT"_recall_reflag.sam

# Extract unmapped reads with samtools
samtools view -bf 4 \
                "$OUTPUT"_recall_reflag.sam \
                > "$OUTPUT"_recall_reflag_unmapped.bam

# Convert bam file to fastq
samtools fastq "$OUTPUT"_recall_reflag_unmapped.bam \
             > "$OUTPUT"_recall_reflag_unmapped.fastq

# Delete intermediate files
rm $OUTPUT.bam
rm $OUTPUT.sam
rm "$OUTPUT"_recall_reflag_unmapped.bam
