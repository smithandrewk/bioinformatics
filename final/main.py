#!/usr/bin/env python3
# author : Andrew Smith
# date : 120820 @ 20:02
# file : main.py
# description : to get proteome from NCBI downloadable CSV and determine which protein any given residue is in

import sys
from Bio import SeqIO

def get_seq(protein):
    """

    """
    for record in SeqIO.parse("SARS_CoV_2_Protein_Sequences.fasta","fasta"):
        if(record.id==protein["Protein product"]):
            return record.seq
    print("not found")

def get_proteome(infile):
    file = open(infile,'r')
    header = [] # ['Name', 'Accession', 'Start', 'Stop', 'Strand', 'GeneID', 'Locus', 'Locus tag', 'Protein product', 'Length', 'Protein Name']
    proteome = []
    for line in file:
        if(line[0]=="#"): # header
            header = line[1:].split(',')
            header[len(header)-1] = header[len(header)-1][:-1]
        else:
            protein = {}
            for i,item in enumerate(line.split(',')):
                protein[header[i]] = item.replace("\"","")
            protein[header[len(header)-1]] = protein[header[len(header)-1]][:-1]
            proteome.append(protein)
    file.close()
    return header,proteome

def which_protein(residue):
    for protein in proteome:
        if(residue<int(protein["Start"])):
            continue
        elif (residue<int(protein["Stop"])):
            return protein
    return "None"

if(len(sys.argv)!=2):
    print("Usage: ./main.py <infile>")
    sys.exit()

infile = sys.argv[1]

header,proteome = get_proteome(infile)
protein = which_protein(500)
count = {}
count["None"]=0
for snp in open("snps","r"):
    if(which_protein(int(snp))=="None"):
        count["None"] += 1
        continue
    if(which_protein(int(snp))["Protein product"] not in count):
        count[which_protein(int(snp))["Protein product"]] = 1
    else:
        count[which_protein(int(snp))["Protein product"]] += 1
for key in count:
    print(key,count[key])
