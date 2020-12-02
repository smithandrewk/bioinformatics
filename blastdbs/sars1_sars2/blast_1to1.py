#!/usr/bin/env python3
# author : Andrew Smith
# date : 111820 @ 12:57
# file : blast.py
# description : provide initial statistics of a query against a blast searchable database
# NOTE: this script is specifically configured for one query and one subject.
import sys # command line arguments
from time import time # execution time
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from numpy import mean

if(len(sys.argv)!=3):
    print("Usage : ./main.py <DBname> <QueryFile>")
    sys.exit()

start_time = time()

dbname = sys.argv[1]
queryfile = sys.argv[2]

blastnCommandLine = NcbiblastnCommandline(query=queryfile,db=dbname,outfmt=5,out="results.xml")
stdout,stderr = blastnCommandLine()

num_queries = 0 # number of sequences in query
num_alignments = 0 # number of sequence in blast searchable database
num_hsps = 0 # number of high scoring pairs in alignment
total_length = 0
total_identities = 0
total_gaps = 0
es = []
for query in NCBIXML.parse(open("results.xml")):
    num_queries += 1
    for alignment in query.alignments:
        num_alignments += 1
        for hsp in alignment.hsps:
            print("___________HSP: %d____________" % num_hsps)
            num_hsps += 1
            l = hsp.align_length # note that len(hsp.sbjct) == len(hsp.query) == hsp.align_length
            total_length += l
            total_identities += hsp.identities
            total_gaps += hsp.gaps
            iden_percent = hsp.identities/l
            gap_percent = hsp.gaps/l
            print("E : "+str(hsp.expect)+"%") # percent chance that this high-scoring segment pair would occur by chance
            print("Identities : %d" % hsp.identities +"/"+str(l)+" => "+str(iden_percent)+"%")
            print("Gaps : %d" % hsp.gaps+"/"+str(l)+" => "+str(gap_percent)+"%")
            es.append(hsp.expect)
identities_weighted = total_identities/total_length
gaps_weighted = total_gaps/total_length
print("___________Summary___________")
print("Expectation: "+str(mean(es)))
print("Identities Weighted : %.2f" % (identities_weighted*100)+"%")
print("Gaps Weighted : %.2f" % (gaps_weighted*100)+"%")
print("Num_HSPs : ",num_hsps)
print("Num_Alignments : ",num_alignments)
print("Num_Queries : ",num_queries)
print("Execution Time : ",time()-start_time,"s")
