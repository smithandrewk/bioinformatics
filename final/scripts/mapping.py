#!/usr/bin/env python3
from Bio import SeqIO
from Bio import AlignIO
from main import *

sars = SeqIO.read("Sars_CoV_2","fasta") # get unaligned sars
sars_align = AlignIO.read("alignment","fasta")._records[0] # get aligned sars
header,proteome = get_proteome("proteome.csv") # get unaligned proteome

mapping = [0]*len(sars_align) # init container for mapping
j = 0
for i,nucl in enumerate(sars_align): # for each nucleotide in aligned sars
    if(nucl=="-"): # if indel
        mapping[i] = -1 # don't map
    else: # if not
        mapping[i] = j # map
        j += 1 # iterate

for protein in proteome: # for each protein in unaligned proteome
    start = int(protein["Start"]) # get start for this protein
    stop = int(protein["Stop"]) # get stop for this protein
    protein["Start"] = mapping.index(start)
    protein["Stop"] = mapping.index(stop)

file = open("unaligned_proteome.csv","w")
file.write("#"+",".join(header)+"\n")
for protein in proteome:
    values = []
    for key in protein:
        values.append(str(protein[key]))
    file.write(",".join(values)+"\n")
file.close()

