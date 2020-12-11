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

def get_protein_seq_from_name(protein,proteome):
    """

    """
    record = SeqIO.read("Sars_CoV_2","fasta")
    sars = record.seq
    start = proteome[protein]["Start"]-1
    stop = proteome[protein]["Stop"]
    # print(start,stop)
    # print(stop-start)
    return sars[start:stop].translate()

def get_proteome(filename):
    file = open(filename,"r")
    proteome = {}
    for line in file:
        if (line[0]=="#"):
            #header
            continue
        else:
            protein = {}
            split = line.replace("\n","").replace("\"","").split(",")
            protein["Start"]= int(split[1])
            protein["Stop"]=int(split[2])
            proteome[split[0]]=protein
    return proteome

def which_protein(residue,proteome):
    for protein in proteome:
        if(residue<proteome[protein]["Start"]):
            continue
        elif (residue<proteome[protein]["Stop"]):
            return protein
    return "none"

def get_protein_count_for_snps(infile,proteome):
    count = {}
    count["None"]=0
    for snp in open(infile,"r"):
        if(which_protein(int(snp),proteome)=="None"):
            count["None"] += 1
            continue
        if(which_protein(int(snp),proteome) not in count):
            count[which_protein(int(snp),proteome)] = 1
        else:
            count[which_protein(int(snp),proteome)] += 1
    return count

def get_codon_wrt_protein(protein,snp,proteome):
    start = proteome[protein]["Start"]
    aa_index = (snp-start)//3
    codon_index = (snp-start) % 3
    nucl_index = (aa_index*3)+start-1 # minus 1 because seq starts at index 0
    record = SeqIO.read("Sars_CoV_2","fasta")
    alignment = AlignIO.read("data/alignment","fasta")
    record = alignment._records[0]
    codon = record.seq[nucl_index:nucl_index+3]
    # print_codon(codon,codon_index)
    return codon,codon_index

def print_codon(codon,codon_index):
    for i,nucl in enumerate(codon):
        if(i==codon_index):
            print(bcolors.FAIL + nucl + bcolors.ENDC,end='')
        else:
            print(nucl,end='')
    print("=>"+codon.translate())
# def get_proteome_from_csv(infile):
#     file = open(infile,'r')
#     header = [] # ['Name', 'Accession', 'Start', 'Stop', 'Strand', 'GeneID', 'Locus', 'Locus tag', 'Protein product', 'Length', 'Protein Name']
#     proteome = []
#     for line in file:
#         if(line[0]=="#"): # header
#             header = line[1:].split(',')
#             header[len(header)-1] = header[len(header)-1][:-1]
#         else:
#             protein = {}
#             for i,item in enumerate(line.split(',')):
#                 protein[header[i]] = item.replace("\"","")
#             protein[header[len(header)-1]] = protein[header[len(header)-1]][:-1]
#             proteome.append(protein)
#     file.close()
#     return header,proteome