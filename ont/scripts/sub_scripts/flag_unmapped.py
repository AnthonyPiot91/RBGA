"""
Reflag mapped reads with specified clipped bases attributes in a SAM file

Usage: python flag_unmapped.py ref.fasta.fai input.sam output.sam
"""

import re
import sys

# Minimum mapping quality
MIN_MAPQUAL = 30
# Maximun percentage of sequence divergence
MAX_SEQ_DIV = 4
# Minimum number of aligned bases required to keep read
MIN_MATCH = 1000
# Maximum number of clipped bases at a read extremity and within the reference
MAX_CLIP_WITHIN = 100
# Maximum number of clipped bases at a read extremity and within the reference
# for reads going past the reference
MAX_CLIP_WITHIN_EDGE = 200

# Define input and output file
REF_FILE = open(sys.argv[1], "r")
SAM_FILE = open(sys.argv[2], "r")
OUT_SAM = open(sys.argv[3], "w")

# Define function to reflag reads as "unmapped"
def reflag():
    global line
    col[1] = "4"
    col[2] = "*"
    col[3] = "0"
    col[4] = "0"
    col[5] = "*"
    if len(col) == 23:
        del col[20] # remove SA flag
    del col[-1] # remove MD flag
    line = "\t".join(col) + "\n"

# Create a dictionary to store the name and length of each
# sequence in the reference
dict = {}
for line_ref in REF_FILE:
    col_ref = line_ref.split("\t")
    seq_ids = col_ref[0]
    length = col_ref[1]
    dict[seq_ids] = length

# Go through each line in the input SAM file
for line in SAM_FILE:

    # Write headers to SAM output
    if line.startswith("@"):
        OUT_SAM.write(line)

    # Go through each read in the SAM input
    else:
        col = line.split("\t")

        # Store read attributes
        flag = col[1]
        id = col[2]
        read_start = int(col[3])
        mapqual = col[4]
        CIGAR_str = col[5]
        MD = col[-1]

	    # If unmapped read write to output and continue to next read
        if flag == "4":
            OUT_SAM.write(line)
            continue

        # Remove secondary alignment
        elif flag == "256" or flag == "272":
            continue

        # Remove supplementary alignment
        elif flag == "2048" or flag == "2064":
            continue

        # Reflag read if mapping quality < to MIN_MAPQUAL
        elif int(mapqual) < MIN_MAPQUAL:
            reflag()

        else:

            # Store length of the reference sequence
            ref_end = int(dict[id])

            # Calculate the number of aligned bases
            match = re.findall("[0-9]*M", CIGAR_str)
            match = [i.replace("M","") for i in match]
            match = [int(i) for i in match]
            sum_match = sum(match)

            # Reflag read if the number of aligned bases is less than MIN_MATCH
            if sum_match < MIN_MATCH:
                reflag()

            # Reflag read if a stretch of more than 99 insertions, deletions or
            # clipped bases inside the read
            elif re.search("[A-Z][0-9]{3,}[I,D,H,S][0-9]", CIGAR_str):
                reflag()

            else:

                # Calculate number of insertion
                matchI = re.findall("[0-9]*I", CIGAR_str)
                matchI = [i.replace("I","") for i in matchI]
                matchI = [int(i) for i in matchI]
                num_Ins = sum(matchI)

                # Calculate number of deletion
                matchD = re.findall("[0-9]*D", CIGAR_str)
                matchD = [i.replace("D","") for i in matchD]
                matchD = [int(i) for i in matchD]
                num_Del = sum(matchD)

                # Calculate number of mismatch
                matchM = re.findall("[0-9]A|[0-9]T|[0-9]G|[0-9]C", MD)
                matchM = [re.sub("[0-9]","", i) for i in matchM]
                num_Mis = len(matchM)

                # Calculate aligned read length
                matchL = re.findall("[0-9]+", MD)
                matchL = [int(i) for i in matchL]
                aligned_read_len = sum(matchL) + num_Ins +num_Mis

                # Reflag read if sequence divergence between read and ref > MAX_SEQ_DIV
                total_mismatch = num_Ins + num_Del + num_Mis
                seq_div = (total_mismatch * 100) / aligned_read_len
                if seq_div > MAX_SEQ_DIV:
                    reflag()

                # If more than 99 clipped bases at the beggining of the read
                elif re.findall("^[0-9]{3,}[H,S]", CIGAR_str):
                    match2 = re.findall("^[0-9]{3,}[H,S]", CIGAR_str)
                    clipped_bases_start = match2[0]
                    clipped_bases_start = clipped_bases_start.replace("H","").replace("S","")
                    clipped_bases_start = int(clipped_bases_start)

                    # Reflag read if clipped bases are within the reference sequence
                    # and greater than MAX_CLIP_WITHIN
                    if clipped_bases_start < read_start and clipped_bases_start > MAX_CLIP_WITHIN:
                        reflag()

                    # Reflag read if clipped bases reach the start of the reference and more than
                    # MAX_CLIP_WITHIN_EDGE clipped bases are within the reference sequence
                    elif clipped_bases_start > read_start and read_start > MAX_CLIP_WITHIN_EDGE:
                        reflag()

                # If more than 99 clipped bases at the end of the read
                elif re.findall("[0-9]{3,}[H,S]$", CIGAR_str):
                    match2 = re.findall("[0-9]{3,}[H,S]$", CIGAR_str)
                    clipped_bases_end = match2[0]
                    clipped_bases_end = clipped_bases_end.replace("H","").replace("S","")
                    clipped_bases_end = int(clipped_bases_end)

                    # Calculate the number of aligned bases and inserted gaps
                    match3 = re.findall("[0-9]*[I,M]", CIGAR_str)
                    match3 = [i.replace("M","").replace("I","") for i in match3]
                    match3 = [int(i) for i in match3]
                    sum_match3 = sum(match3)
                    aligned_read_end = read_start + sum_match3
                    clipped_end_pos = aligned_read_end + clipped_bases_end

                    # Reflag read if clipped bases are within the reference sequence and greater
                    # than MAX_CLIP_WITHIN
                    if clipped_end_pos < ref_end and clipped_bases_end > MAX_CLIP_WITHIN:
                        reflag()

                    # Reflag read if clipped bases reach the end of the reference and more than
                    # MAX_CLIP_WITHIN_EDGE clipped bases are within the reference sequence
                    else:
                        clipped_end_within = ref_end - aligned_read_end
                        if clipped_end_pos > ref_end and clipped_end_within > MAX_CLIP_WITHIN_EDGE:
                            reflag()

            OUT_SAM.write(line)

OUT_SAM.close()
