#!/usr/bin/env python3
# author : Andrew Smith
# date : 120820 @ 20:02
# file : main.py
# description : to get proteome from NCBI downloadable CSV and determine which protein any given residue is in

import sys
from Bio import SeqIO
from Bio import AlignIO

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
    for record in SeqIO.parse("data/sequences.fasta","fasta"):
        if(record.id==protein["Protein product"]):
            return record.seq
    print("not found")

def get_proteome_from_csv(infile):
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

def which_protein(residue,proteome):
    for protein in proteome:
        if(residue<int(protein["Start"])):
            continue
        elif (residue<int(protein["Stop"])):
            return protein
    return "None"

def get_protein_for_each_snp(infile,proteome):
    count = {}
    count["None"]=0
    for snp in open(infile,"r"):
        if(which_protein(int(snp),proteome)=="None"):
            count["None"] += 1
            continue
        if(which_protein(int(snp),proteome)["Protein product"] not in count):
            count[which_protein(int(snp),proteome)["Protein product"]] = 1
        else:
            count[which_protein(int(snp),proteome)["Protein product"]] += 1
    for key in count:
        print(key,count[key])

def get_codon_wrt_protein(protein,snp):
    aa_index = (snp-int(protein["Start"]))//3
    codon_index = (snp-int(protein["Start"])) % 3
    nucl_index = (aa_index*3)+int(protein["Start"])
    #record = SeqIO.read("data/Sars_CoV_2","fasta")
    alignment = AlignIO.read("data/alignment","fasta")
    record = alignment._records[0]
    codon = record.seq[nucl_index:nucl_index+3]
    print(aa_index)
    #aa = get_seq(protein)[aa_index]
    return codon, codon_index
def wtf(proteome):
    Sars_CoV_2 = SeqIO.read("data/Sars_CoV_2","fasta")
    for snp in open("data/snps","r"):    
        print("____________________")
        snp = int(snp)
        print(snp)
        protein = which_protein(snp,proteome)
        if(protein == "None"):
            print("Not in protein")
        else:
            codon, codon_index = get_codon_wrt_protein(protein,snp)
            for i,nucl in enumerate(codon):
                if(i==codon_index):
                    print(bcolors.FAIL + nucl + bcolors.ENDC,end='')
                else:
                    print(nucl,end='')
            print()
