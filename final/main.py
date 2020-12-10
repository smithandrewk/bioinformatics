#!/usr/bin/env python3
# author : Andrew Smith
# date : 120820 @ 20:02
# file : main.py
# description : to get proteome from NCBI downloadable CSV and determine which protein any given residue is in

import sys
from Bio import SeqIO

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
def get_seq(protein):
    """

    """
    for record in SeqIO.parse("sequences.fasta","fasta"):
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
def get_protein_for_each_snp():
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
def get_codon_wrt_protein(protein,snp):
    aa_index = (snp-int(protein["Start"]))//3
    codon_index = (snp-int(protein["Start"])) % 3
    nucl_index = (aa_index*3)+int(protein["Start"])
    record = SeqIO.read("Sars_CoV_2","fasta")
    codon = record.seq[nucl_index:nucl_index+3]
    print(aa_index)
    aa = get_seq(protein)[aa_index-1]
    return codon, codon_index,aa

#if(len(sys.argv)!=3):
#    print("Usage: ./main.py <infile> <snp index>")
#    sys.exit()
#
#infile = sys.argv[1]
#snp = int(sys.argv[2])
#
#header,proteome = get_proteome(infile)
#Sars_CoV_2 = SeqIO.read("Sars_CoV_2","fasta")
#for snp in open("snps","r"):    
#    print("____________________")
#    snp = int(snp)
#    print(snp)
#    protein = which_protein(snp)
#    if(protein == "None"):
#        print("Not in protein")
#    else:
#        codon, codon_index,aa = get_codon_wrt_protein(protein,snp)
#        for i,nucl in enumerate(codon):
#            if(i==codon_index):
#                print(bcolors.FAIL + nucl + bcolors.ENDC,end='')
#            else:
#                print(nucl,end='')
#        print("-->",aa)
