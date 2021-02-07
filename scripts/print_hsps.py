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

if(len(sys.argv)!=3):
    print("Usage : ./main.py <DBname> <QueryFile>")
    sys.exit()

start_time = time()

dbname = sys.argv[1]
queryfile = sys.argv[2]

blastnCommandLine = NcbiblastnCommandline(query=queryfile,db=dbname,outfmt=5,out="results.xml")
stdout,stderr = blastnCommandLine()

for query in NCBIXML.parse(open("results.xml")):
    for alignment in query.alignments:
        for hsp in alignment.hsps:
            print(hsp)
