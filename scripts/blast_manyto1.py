#!/usr/bin/env python3
# author : Andrew Smith
# date : 112020 @ 18:09
# file : blast_manyto1.py
# description : 
import sys # command line arguments
from time import time # execution time
from Bio.Blast.Applications import NcbiblastnCommandline # to format blast command line
from Bio.Blast import NCBIXML # to parse blast output xml file
import matplotlib.pyplot as plt
from numpy import mean

if(len(sys.argv)!=3):
    print("Usage : ./main.py <DBname> <QueryFile>")
    sys.exit()

start_time = time()

dbname = sys.argv[1]
queryfile = sys.argv[2]

blastnCommandLine = NcbiblastnCommandline(query=queryfile,db=dbname,outfmt=5,out="results.xml")
stdout,stderr = blastnCommandLine()

mutations = []
hsps = []
es = []
for query in NCBIXML.parse(open("results.xml")): # for each query sequence
    for alignment in query.alignments: # for every sequence in the blast searchable database
        hsps.append(len(alignment.hsps)) 
        for hsp in alignment.hsps: # for every high-scoring segment pair
            query_start = hsp.query_start # get starting residue of query sequence
            index = 23403-query_start # find adjusted mutation index of desired mutation
            if(index <0 or index > len(hsp.query)): # mutation index not contained in this hsp
                continue
            else: # mutation contained in this hsp
                mutations.append(hsp.sbjct[index]) # append to mutations
                es.append(hsp.expect)
                # TODO : this magically produces one mutation per alignment. HSPs, ideally, agree on the mutation; therefore, we should check whether they agree, then append for each alignment
for nucl in set(mutations): # for each type of nucleotide in mutations
    print("%.2f" % (mutations.count(nucl)/len(mutations)*100)+"%",nucl) # print percentage of this nucleotide
print(mean(hsp.expect))

print("\nExecution Time : ",time()-start_time)

plt.hist(mutations)
plt.title("A23403X : X in 193 Sars-CoV-2 Sequences from South Carolina")
plt.xlabel("X aligned to residue 23403 in Sars-CoV-2 Reference Sequence")
plt.ylabel("Frequency")
plt.savefig("a23403x.jpg")
plt.show()
