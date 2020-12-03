#!/usr/bin/env python3
# author : Andrew Smith
# date : 120220 @ 00:44
# file : analyze.py
# description : base script along with plotting fig 1
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
for i in tqdm(range(len(parent.seq))):
    num_muts.append(0)
    for record in alignment._records[1:]:
        if(parent.seq[i]!=record.seq[i]):
            num_muts[i]+=1
    perc_muts.append(num_muts[i]/len(alignment._records[1:])*100)
snps = [i for i,perc_mut in enumerate(perc_muts) if perc_mut>5]
## Fig 1
plt.ylabel("Percentage Mutation")
plt.xlabel("Residue")
plt.axhline(5,color="r",linestyle="--")
plt.plot(perc_muts)
plt.savefig("fig1.jpg",dpi=100,bbox_inches='tight')
plt.show()
