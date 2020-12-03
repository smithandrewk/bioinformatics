#!/usr/bin/env python3
# author : Andrew Smith
# date : 120220 @ 00:42
# file : clustalo.py
# description : a useless python script which runs ClustalOmega
import sys
from Bio.Align.Applications import ClustalOmegaCommandline

if(len(sys.argv)!=3):
    print("Usage: ./main.py <infile> <outfile>")
    sys.exit()

cmd = ClustalOmegaCommandline("clustalo",infile=sys.argv[1],outfile=sys.argv[2])

stdout,stderr = cmd()



