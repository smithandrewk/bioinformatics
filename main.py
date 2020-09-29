#!/usr/bin/env python3
#author: Andrew Smith
#last edited 092420 @ 10:37
#description: to read in a fasta file

import sys
import time
from structures import Sequence, Gene
from read_fasta import *

#Proper Script Usage
if(len(sys.argv)!=4):
    print("Usage: ./read_fasta.py {infile} {outfile} {searchpattern}")
    sys.exit(0)
#Parse Command Line Arguments
infile = str(sys.argv[1])
outfile = str(sys.argv[2])
searchpattern = str(sys.argv[3])

start_time = time.time()# start execution timer

# We obtain parent sequence, write it to file, and keep in memory
# because we compare each sequence against it
parent = get_next_sequence(infile,outfile,"w+",searchpattern)
# scan_for_genes(parent)
last_line = parent.last_line

while(last_line!=-1): # while file is not empty
    print(last_line)
    seq = get_next_sequence(infile,outfile,"a",searchpattern,last_line,parent) # read next sequence
    last_line = seq.last_line

print("--- %s seconds ---" % (time.time() - start_time))


