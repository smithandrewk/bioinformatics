#!/usr/bin/env python3
import sys
from time import time
from tqdm import tqdm
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


start_time = time()
# get cov1
infile = open("seq1", "r")
seqHandle = SeqIO.parse(infile, "fasta")
cov1 = next(seqHandle).seq
infile.close()
# get cov2
infile = open("seq2", "r")
seqHandle = SeqIO.parse(infile, "fasta")
cov2 = next(seqHandle).seq
infile.close()
# open outfile
outfile = open("out", "w")
# write description
outfile.write(">Optimal Global Alignments for SarsCov1 and SarsCov2 with match = 1, mistmatch = -1, gap = -8:\n")
# perform alignment
alignment = pairwise2.align.globalms(cov1, cov2, 1, -1, -8, -8)
print("Alignment Time : ",time()-start_time)
for align in tqdm(alignment):
    outfile.write(format_alignment(*align))
outfile.close()
print("Execution Time: ", time() - start_time)
