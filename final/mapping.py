#!/usr/bin/env python3
from Bio import SeqIO
from Bio import AlignIO
from main import *
sars = SeqIO.read("Sars_CoV_2","fasta")
sars_align = AlignIO.read("alignment","fasta")._records[0]
header,proteome = get_proteome("proteome.csv")
gap = sars_align.seq.find("-")
mapping = [0]*len(sars_align)
count = 0
for i,nucl in enumerate(sars_align):
    if(nucl=="-"):
        mapping[i] = -1
    else:
        mapping[i]=count
        count += 1
for protein in proteome:
    start = int(protein["Start"])
    stop = int(protein["Stop"])
    print("__________________")
    print(protein["Start"],protein["Stop"])
    print(mapping.index(start),mapping.index(stop))


