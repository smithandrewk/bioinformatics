import sys
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices as matlist
from time import time

if __name__ == '__main__':
    start_time = time()
    infile = open(sys.argv[1], "r")
    outfile = open(sys.argv[2],"w")
    seqHandle = SeqIO.parse(infile,"fasta")
    matrix = matlist.load("BLOSUM50")
    S0 = next(seqHandle)
    S1 = next(seqHandle)
    S2 = next(seqHandle)
    outfile.write(">Opimal Global Alignments for S0 and S1 with match = 1, mismatch = -1, gap = -8:\n")
    alignment = pairwise2.align.globalms(S0.seq, S1.seq,1,-1,-8,-8)
    for align in alignment:
         outfile.write(format_alignment(*align))
    outfile.write(">Optimal Global Alignment for S0 and S2 with match = 1, mismatch = -1, gap = -8:\n")
    alignment = pairwise2.align.globalms(S0.seq, S2.seq,1,-1,-8,-8)
    for align in alignment:
         outfile.write(format_alignment(*align))
    print("Execution Time: ",time()-start_time)
    infile.close()
    outfile.close()
