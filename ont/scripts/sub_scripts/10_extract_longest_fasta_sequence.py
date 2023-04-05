"""
Sort a fasta file by decreasing sequence length and extract the first (longest) sequence.

Used to extract the chromosome contig from ONT assemblies of Borrelia strains.

Usage: python this_script fasta_file_to_sort longest_sequence_output
"""

import sys
from Bio import SeqIO

record = list(SeqIO.parse(sys.argv[1], "fasta"))

record = sorted(record, key = len, reverse = True)

if len(record[0]) < 850000:
	print("Warning, your longest contig is smaller than 850'000bp!")

SeqIO.write(record[0], sys.argv[2], "fasta")
