# Define parameters
MAPPING=$1
REF=$2
FASTQ=$3
OUTPUT=$4
THREADS=$5
RACON_ERR_THR=$6
SAM=$7

# Print racon version
racon --version

# Run racon
racon   $FASTQ \
        $SAM \
        $REF \
        -t $THREADS \
        --error-threshold $RACON_ERR_THR \
        --no-trimming \
        > "$OUTPUT"_racon.fa
