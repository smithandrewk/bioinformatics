#!/usr/bin/env python3
# author : Andrew Smith
# date : 120220 @ 00:53
# file : analyze.py
# description : fig 3 script but this produces the snps table file
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

## SNP Table
breakdown = [] # container for breakdown of each snp
for snp_mut in snp_muts: # for each snp
    ap = { # "ap" meaning "to append"
        "A":0,
        "T":0,
        "C":0,
        "G":0,
        "-":0,
        "N":0
        }
    for mut in snp_mut: # for each mutation in this snp
        ap[mut]+=1 # iterate the count of that mutation
    for i in ap: # then correct for number of sequences and make percentage
        ap[i]=ap[i]/193*100
    breakdown.append(ap) # append to list of breakdowns
# Writing to file
outfile = open("snps","w")
for i,snp in enumerate(snps):
    outfile.write(str(snp)+","+str(perc_muts[snp])+","+str(breakdown[i])+"\n")
outfile.close()

