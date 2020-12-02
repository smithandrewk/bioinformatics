#!/usr/bin/env python3

import sys
from Bio.Align.Applications import ClustalOmegaCommandline

if(len(sys.argv)!=3):
    print("Usage: ./main.py <infile> <outfile>")
    sys.exit()

cmd = ClustalOmegaCommandline("clustalo",infile=sys.argv[1],outfile=sys.argv[2])

stdout,stderr = cmd()



