#!/usr/bin/env python3
#to test the randomness, or sudorandomness, of the choice function from the random package in python in the context of bioinformatics
from random import choice
nucleic_acids = ['A','T','C','G']
a=0
t=0
c=0
g=0
o=0
for i in range(10000):
    ch = choice(nucleic_acids)
    if(ch=='A'):
        a=a+1
    elif(ch=='T'):
        t=t+1
    elif(ch=='C'):
        c=c+1
    elif(ch=='G'):
        g=g+1
    else:
        o=o+1
print(a,t,c,g,o)
