# Define parameters
MAPPING=$1
REF=$2
FASTQ=$3
OUTPUT=$4

# Print minimap2 version
minimap2 --version

# Run minimap2
minimap2 -a \
         -x $MAPPING \
         $REF \
         $FASTQ \
         -o $OUTPUT.sam
