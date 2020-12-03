#!/usr/bin/env python3
# author : Andrew Smith
# date : 120220 @ 00:46
# file : analyze.py
# description : different from the base script in that we record the actual mutations for each sequence in the SNPs, then plotting fig3
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

num_muts = []
perc_muts = []
snp_muts = []
snps = []
for i in tqdm(range(len(parent.seq))):
    num_muts.append(0)
    muts = []
    for record in alignment._records[1:]:
        muts.append(record.seq[i])
        if(parent.seq[i]!=record.seq[i]):
            num_muts[i]+=1
    perc_mut = num_muts[i]/len(alignment._records[1:])*100
    perc_muts.append(perc_mut)
    if(perc_mut>5):
        snp_muts.append(muts)
        snps.append(i)
## Fig 3
content = {
        "A":0,
        "T":0,
        "C":0,
        "G":0,
        "-":0,
        "N":0
        }
for snp_mut in snp_muts:
    for mut in snp_mut:
        content[mut]+=1
freqs = [content[i]/(len(snp_muts)*193) for i in content]
plt.bar(content.keys(),freqs)
plt.ylabel("Frequency")
plt.savefig("fig3.jpg",dpi=100,bbox_inches='tight')
plt.show()
