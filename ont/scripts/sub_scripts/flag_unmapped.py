"""
Reflag mapped reads with specified clipped bases attributes in a SAM file

Usage: python flag_unmapped.py ref.fasta.fai input.sam output.sam
"""

import re
import sys

# Minimal number of aligned bases required to keep read
MATCH = 1000
# Maximum number of clipped bases at a read extremity and within the reference
MAX_CLIP_WITHIN = 200
# Maximum number of clipped bases at a read extremity and within the reference
# for reads going past the reference
MAX_CLIP_WITHIN_EDGE = 200

REF_FILE = open(sys.argv[1], "r")
SAM_FILE = open(sys.argv[2], "r")
OUT_SAM = open(sys.argv[3], "w")

# Create a dictionary to store the name and length of each
# sequence in the reference
dict = {}
for line_ref in REF_FILE:
    col_ref = line_ref.split("\t")
    seq_ids = col_ref[0]
    length = col_ref[1]
    dict[seq_ids] = length

# Go through each line in input SAM file
for line in SAM_FILE:

    # Write headers to SAM output
    if line.startswith("@"):
        OUT_SAM.write(line)

    # Go through each aligned reads in the SAM input
    else:
        col = line.split("\t")

        # If unmapped read write to output and continue to next read
        if col[1] == "4":
            OUT_SAM.write(line)
            continue

        # Remove secondary alignment
        elif col[1] == "256" or col[1] == "272":
            continue

        # Remove supplementary alignment
        elif col[1] == "2048" or col[1] == "2064":
            continue

        # Store read attributes
        id = col[2]
        read_start = int(col[3])
        CIGAR_str = col[5]

        # Store reference sequence attribute
        ref_end = int(dict[id])

        # Calculate the number of aligned bases
        match = re.findall("[0-9]*M", CIGAR_str)
        match = [i.replace("M","") for i in match]
        match = [int(i) for i in match]
        sum_match = sum(match)

        # Reflag if the number of aligned bases is less than MATCH
        if sum_match < MATCH:
            col[1] = "4"
            col[2] = "*"
            col[3] = "0"
            col[4] = "0"
            col[5] = "*"
            line = "\t".join(col)

        # Reflag if a stretch of more than 99 insertions or deletions inside the read
        elif re.search("[A-Z][0-9]{3,}[D,I][0-9]", CIGAR_str):
            col[1] = "4"
            col[2] = "*"
            col[3] = "0"
            col[4] = "0"
            col[5] = "*"
            line = "\t".join(col)

        # Reflag if a strech of more than 99 clipped bases inside the read
        elif re.search("[A-Z][0-9]{3,}[H,S][0-9]", CIGAR_str):
            col[1] = "4"
            col[2] = "*"
            col[3] = "0"
            col[4] = "0"
            col[5] = "*"
            line = "\t".join(col)

        else:

            # If more than 99 clipped bases at the beggining of the read
            if re.findall("^[0-9]{3,}[H,S]", CIGAR_str):
                match2 = re.findall("^[0-9]{3,}[H,S]", CIGAR_str)
                clipped_bases_start = match2[0]
                clipped_bases_start = clipped_bases_start.replace("H","").replace("S","")
                clipped_bases_start = int(clipped_bases_start)

                # Reflag if clipped bases are within the reference sequence
                # and greater than MAX_CLIP_WITHIN
                if clipped_bases_start < read_start and clipped_bases_start > MAX_CLIP_WITHIN:
                    col[1] = "4"
                    col[2] = "*"
                    col[3] = "0"
                    col[4] = "0"
                    col[5] = "*"
                    line = "\t".join(col)

                # Reflag if clipped bases reach the start of the reference and more than
                # MAX_CLIP_WITHIN_EDGE clipped bases are within the reference sequence
                elif clipped_bases_start > read_start and read_start > MAX_CLIP_WITHIN_EDGE:
                    col[1] = "4"
                    col[2] = "*"
                    col[3] = "0"
                    col[4] = "0"
                    col[5] = "*"
                    line = "\t".join(col)

            # If more than 99 clipped bases at the end of the read
            if re.findall("[0-9]{3,}[H,S]$", CIGAR_str):
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

                # Reflag if clipped bases are within the reference sequence and greater
                # than MAX_CLIP_WITHIN
                if clipped_end_pos < ref_end and clipped_bases_end > MAX_CLIP_WITHIN:
                    col[1] = "4"
                    col[2] = "*"
                    col[3] = "0"
                    col[4] = "0"
                    col[5] = "*"
                    line = "\t".join(col)

                # Reflag if clipped bases reach the end of the reference and more than
                # MAX_CLIP_WITHIN_EDGE clipped bases are within the reference sequence
                else:
                    clipped_end_within = ref_end - aligned_read_end
                    if clipped_end_pos > ref_end and clipped_end_within > MAX_CLIP_WITHIN_EDGE:
                        col[1] = "4"
                        col[2] = "*"
                        col[3] = "0"
                        col[4] = "0"
                        col[5] = "*"
                        line = "\t".join(col)
        OUT_SAM.write(line)

OUT_SAM.close()
