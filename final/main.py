#!/usr/bin/env python3
from utils import *
import operator
from Bio.Seq import Seq

if(len(sys.argv)!=1):
    print("Usage: ./main.py")
    sys.exit()


# SNPS START FROM INDEX 1
proteome = get_proteome("aligned_proteome.csv")
# count = get_protein_count_for_snps("snps",proteome)


for lin in open("summary","r"):
    ## parse summary file
    spl = lin.index(",")
    snp = int(lin[0:spl])

    lin = lin[spl+1:]
    spl = lin.index(",")
    mutation_rate = float(lin[0:spl])

    breakdown = eval(lin[spl+1:])

    # Which protein is it in?
    protein = which_protein(snp,proteome)
    # If SNP not in a protein
    if(protein == "none"):
        print("Not in protein")
    else:
        # Get the protein sequence
        seq = get_protein_seq_from_name(protein,proteome)
        # Otherwise, get the codon with respect to the protein
        codon,codon_index = get_codon_wrt_protein(protein,snp,proteome)
        str_codon = str(codon)
        s = list(str_codon)
        key = max(breakdown.items(), key=operator.itemgetter(1))[0]
        if(key==codon[codon_index]):
            continue
        elif(key=="-"):
            continue
        elif(breakdown[key]==0):
            continue
        else: 
            print("Non-synonymous SNP at "+str(snp))
            s[codon_index]=key
            new_codon = Seq("".join(s))
            if(codon!=new_codon):
                print_codon(codon,codon_index)
                print_codon(new_codon,codon_index)
        # for key in breakdown:
        #     if(key==codon[codon_index]):
        #         continue
        #     elif(key=="-"):
        #         continue
        #     elif(breakdown[key]==0):
        #         continue
        #     else: 
        #         print("Non-synonymous SNP at "+str(snp))
        #         s[codon_index]=key
        #         new_codon = Seq("".join(s))
        #         if(codon!=new_codon):
        #             print_codon(codon,codon_index)
        #             print_codon(new_codon,codon_index)