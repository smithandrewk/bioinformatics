#!/usr/bin/env python3
from utils import *

if(len(sys.argv)!=1):
    print("Usage: ./main.py")
    sys.exit()

file = open("summary","r")
outfile = open("aasdgasdg","w")
for line in file:
    print(str(int(line[:line.index(',')])+1))
    print(line[line.index(',')+1:]+"\n")
    outfile.write(str(int(line[:line.index(',')])+1)+line[line.index(','):])

file.close()
outfile.close()

