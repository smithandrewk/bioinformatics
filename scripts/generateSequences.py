#!/usr/bin/env python3
#author: Andrew Smith
#last edited: 091320 @ 14:52
#description: to produce a parent sequence and subsequent sequences with mutations in fasta format for CSCE 555 at UofSC

import sys
from random import choice
from random import sample
from math import ceil

#-----Proper Script Usage-----
if(len(sys.argv)!=5):
    print("Usage: ./generateSequences.py {length of parent sequence} {number of derivative sequences} {percentage of mutations} {output file name}")
    sys.exit(0)
#-----Init Stack and Argument Checks-----
#length of the parent sequence
l = int(sys.argv[1])
if(l<=0):
    print("Length must be greater than 0")
    sys.exit(0)
#number of derivative sequences to be generated from the parent sequence
m = int(sys.argv[2])
if(m<0):
    print("Cannot produce negative number of derivative sequences")
    sys.exit(0)
#percentage of mutations each child sequence should have when compared to the parent (rounded down)
p = float(sys.argv[3])
if(p<0 or p>100):
    print("Percentage must be in range [0,100]")
    sys.exit(0)
#name of the fasta output file
file_name = sys.argv[4]
#letters to choose from for sequence
nucleic_acids = ['A','U','C','G']
#open file object to write to
file_object = open("../data/"+file_name,"w")
#initialize empty parent sequence
parent_sequence = ''

#-----Parent Sequence-----
#file formatting
file_object.write(">S0\n")
#write parent sequence of length l randomly choosing nucleic acids
for i in range(l):
    parent_sequence=parent_sequence+choice(nucleic_acids)
#split parent sequence every 80 characters and store to facilitate 80-long fasta lines
parent_sequence_split = [parent_sequence[i:i+80] for i in range(0,len(parent_sequence),80)]
#write parent sequence to file
file_object.write("\n".join(parent_sequence_split)+"\n")

#-----Derivative Sequences-----
#number of mutations is the ceiling of the percentage times the length
#NOTE: ceil is used to give a reasonable interpretation of p=(1/k)% mutation in a sequence of length less than k*100 where k is a natural number. This implies that any non zero value of p will always produce at least 1 mutation.
number_of_mutations = ceil((p/100)*l)
#array where the value at index i is equal to i to facilitate picking random indices for mutation
index_list = [i for i in range(l)]
#Generate m derivatives of the parent sequence
for i in range(m):
    #set curr sequence to parent sequence as a list to allow to for mutating specific indices of a string
    curr_sequence_list = list(parent_sequence)
    #sample returns a subset of a certain size, so we obtain the indices to mutate
    mutation_indices = sample(index_list,number_of_mutations)
    #sort for printing of mutated indices in fasta comment line
    mutation_indices.sort()
    #fasta comment formatting
    file_object.write(">S"+str(i+1))
    #for each mutation, write the mutation in fasta comment separated by comma
    #to start comma after first mutation index
    first = True
    for mutation in mutation_indices:
        #if first, do not print comma before mutation index
        if first:
            file_object.write(" "+str(mutation))
            first = False
        #otherwise, print comma before mutation index
        else:
            file_object.write(","+str(mutation))
    #newline after gene comment line
    file_object.write("\n")
    #produce mutation for each index in mutation indices  
    for index in mutation_indices:
        #get the current nucleic acid at this index in the parent sequence
        curr_acid = parent_sequence[index]
        #get the index of the current nucleic acid in the nucleic acids array
        curr_acid_index = nucleic_acids.index(curr_acid)
        #get nucleic acids without the current acid
        remaining_acids = nucleic_acids[:curr_acid_index]+nucleic_acids[curr_acid_index+1:]
        #randomly choose an acid for mutation from remaining acids
        curr_sequence_list[index]=choice(remaining_acids)
    #current sequence as string
    curr_sequence="".join(curr_sequence_list)
    #split current sequence every 80 characters to facilitate 80-long fasta lines
    curr_sequence_split = [curr_sequence[i:i+80] for i in range(0,len(curr_sequence),80)]
    #write current sequence to file
    file_object.write("\n".join(curr_sequence_split)+"\n")
#close file object
file_object.close
