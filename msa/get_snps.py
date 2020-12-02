#!/usr/bin/env python3
import sys
from tqdm import tqdm
from Bio import AlignIO
import matplotlib.pyplot as plt
if(len(sys.argv)!=2):
    print("Usage : ./analyze.py <infile>")
    sys.exit()

infile = sys.argv[1]
alignment = AlignIO.read(infile,"fasta")
parent = alignment._records[0]
percentage_mutations = []
snps = []
for nt_index in tqdm(range(0,len(parent))):
    nts = [record.seq[nt_index] for record in alignment._records[1:]]
    if(len(set(nts))==1):
        continue
    nt_counts = {}
    for nt in nts:
        if nt not in nt_counts:
            nt_counts[nt]=1
        else:
            nt_counts[nt]+=1
    snps.append(
            {'index':nt_index,
                'nt_counts':nt_counts})
for snp in snps:
    matches = snp["nt_counts"][parent.seq[snp["index"]]]/193
    if(matches < .95):
        percentage_mutations.append(1-matches)
plt.plot(percentage_mutations)
plt.show()
        

