#!/usr/bin/env python3
# author : Andrew Smith
# date : 120320 @ 00:41
# file : analyze.py
# description : base script which assumes multiple sequence alignment has already been computed, then gets SNPs with a mutation rate greater than 5%
import sys
from tqdm import tqdm
from Bio import AlignIO
if(len(sys.argv)!=2):
    print("Usage : ./analyze.py <infile>")
    sys.exit()

infile = sys.argv[1]
alignment = AlignIO.read(infile,"fasta")
ref = alignment._records[0]

num_muts = []
perc_muts = []

for i in tqdm(range(len(ref.seq))): # for every residue in ref seq
    num_muts.append(0) # append container for number of mutations
    for record in alignment._records[1:]: # for each record which is not ref seq
        if(ref.seq[i]!=record.seq[i]): # if not equal to ref
            num_muts[i]+=1 # is mutation
    perc_muts.append(num_muts[i]/len(alignment._records[1:])*100) # append percent mutations for every residue
snps = [i for i,perc_mut in enumerate(perc_muts) if perc_mut>5] # if percent mutation greater than 5, append residue

