#!/usr/bin/env python3
from utils import *

if(len(sys.argv)!=1):
    print("Usage: ./main.py")
    sys.exit()

header,proteome = get_proteome_from_csv("data/aligned_proteome.csv")
get_protein_for_each_snp("data/snps",proteome)
wtf(proteome)
