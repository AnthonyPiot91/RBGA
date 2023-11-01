# Define paramaters
THREADS=12
echo $THREADS
STRAIN=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../barcode_correspondence.txt | cut -f1)
echo $STRAIN
PARAM="2000_15"
echo $PARAM
REF="09_flye/"$PARAM"_meta_nano-hq/$STRAIN/assembly.fasta"
echo $REF
INPUT_DIR="06_chopper"
echo $INPUT_DIR
OUTPUT="10_tadpole/$STRAIN"
mkdir -p $OUTPUT

# Load conda environment
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate bbmap

tadpole.sh --version

# Run 5 round of contig extension with tadpole
tadpole.sh	threads=$THREADS \
		extendleft=1000 \
		extendright=1000 \
		trimcircular=t \
		overwrite=true \
		in=$REF \
		extra=$INPUT_DIR/"$STRAIN"_trimmed_"$PARAM".fastq \
		out=$OUTPUT_DIR/"$STRAIN"_tadpole1.fa \
		mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=$OUTPUT_DIR/"$STRAIN"_tadpole1.fa \
                extra=$INPUT_DIR/"$STRAIN"_trimmed_"$PARAM".fastq \
                out=$OUTPUT_DIR/"$STRAIN"_tadpole2.fa \
                mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=$OUTPUT_DIR/"$STRAIN"_tadpole2.fa \
                extra=$INPUT_DIR/"$STRAIN"_trimmed_"$PARAM".fastq \
                out=$OUTPUT_DIR/"$STRAIN"_tadpole3.fa \
                mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=$OUTPUT_DIR/"$STRAIN"_tadpole3.fa \
                extra=$INPUT_DIR/"$STRAIN"_trimmed_"$PARAM".fastq \
                out=$OUTPUT_DIR/"$STRAIN"_tadpole4.fa \
                mode=extend

tadpole.sh      threads=$THREADS \
                extendleft=1000 \
                extendright=1000 \
                trimcircular=t \
                overwrite=true \
                in=$OUTPUT_DIR/"$STRAIN"_tadpole4.fa \
                extra=$INPUT_DIR/"$STRAIN"_trimmed_"$PARAM".fastq \
                out=$OUTPUT_DIR/"$STRAIN"_tadpole5.fa \
                mode=extend
