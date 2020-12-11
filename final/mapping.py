#!/usr/bin/env python3
from Bio import SeqIO
from Bio import AlignIO
from utils import *

sars = SeqIO.read("Sars_CoV_2","fasta") # get unaligned sars
sars_align = AlignIO.read("data/alignment","fasta")._records[0] # get aligned sars
proteome = get_proteome("proteome.csv") # get unaligned proteome

mapping = [0]*len(sars_align) # init container for mapping
j = 1
for i,nucl in enumerate(sars_align): # for each nucleotide in aligned sars
    if(nucl=="-"): # if indel
        mapping[i] = -1 # don't map
    else: # if not
        mapping[i] = j # map
        j += 1 # iterate
for protein in proteome: # for each protein in unaligned proteome
    start = proteome[protein]["Start"] # get start for this protein
    stop = proteome[protein]["Stop"] # get stop for this protein
    proteome[protein]["Start"] = mapping.index(start)
    proteome[protein]["Stop"] = mapping.index(stop)

file = open("askd.csv","w")
file.write("#Name,Start,Stop")
for protein in proteome:
    file.write(protein+","+str(proteome[protein]["Start"]+1)+","+str(proteome[protein]["Stop"]+1)+"\n")

file.close()

