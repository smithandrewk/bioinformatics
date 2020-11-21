#!/usr/bin/env python3
# author: Andrew Smith
# date: 100420 @ 17:44
# file: main.py
# description: to analyze an RNA sequence in fasta format

import sys # for command line arguments
sys.path.append('..')
from time import time # to time program execution
from utils.fasta_utilities import *

#Proper Script Usage
if(len(sys.argv)!=4):
    print("Usage: ./read_fasta.py {infile} {outfile} {searchpattern}")
    sys.exit(0)
#Parse Command Line Arguments
infile = str(sys.argv[1])
outfile = str(sys.argv[2])
searchpattern = str(sys.argv[3])

start_time = time()# start execution timer
# We obtain parent sequence, write it to file, and keep in memory
# because we compare each sequence against it
parent = get_next_sequence(infile,outfile,"w+",searchpattern)
offset = parent.offset

while(offset!=-1): # while file is not empty
    seq = get_next_sequence(infile,outfile,"a",searchpattern,offset,parent) # read next sequence
    offset = seq.offset #update offset to utilize seek method for python fileio

print("--- %s seconds ---" % (time() - start_time))
