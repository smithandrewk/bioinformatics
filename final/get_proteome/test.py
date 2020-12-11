#!/usr/bin/env python3
from utils import *

if(len(sys.argv)!=1):
    print("Usage: ./main.py")
    sys.exit()

file = open("GenBank","r")
features = []
feature = []
for line in file:
    if(line[0]!=" "): #store new feature
        features.append(feature)
        feature = []
        b = line.split()
        feature.append(b[0].replace("\n","").replace("\"",""))
        c = b[1].split("..")

        for d in c:
            feature.append(d.replace("\n","").replace("\"",""))
    else:
        a = line.replace(" ","")
        if(a[0:8]=="/product"):
            c = a[1:].split("=")
            feature[0]=c[1].replace("\n","").replace("\"","")

outfile = open("p","w")

for feature in features:
    outfile.write(",".join(feature)+"\n")
outfile.close()
