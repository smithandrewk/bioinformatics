#!/usr/bin/env python3
#author: Andrew Smith
#last edited 092420 @ 10:37
#description: to read in a fasta file

import sys
import time

#Proper Script Usage
if(len(sys.argv)!=4):
    print("Usage: ./read_fasta.py {infile} {outfile} {searchpattern}")
    sys.exit(0)
infile = str(sys.argv[1])
outfile = str(sys.argv[2])
searchpattern = str(sys.argv[3])
class Sequence:
    def __init__(self):
        self.id = "NONE"
        self.mutation_indices = 0
        self.length = 0
        self.body = ""
        self.last_line = 0
        self.length = 0
        self.recovery = False
        self.query_matches = 0
    def __str__(self):
        if(len(self.mutation_indices)==0):
            return self.id+", "+str(self.length)+", "+str(self.recovery)+", "+str(self.query_matches)
        return self.id+", "+str(self.length)+", "+str(self.recovery)+", "+str(self.query_matches)
class Gene:
    def __init__(self):
        self.start_index = 0
        self.stop_index = 0
    def __str__(self):
        return "("+str(self.start_index)+","+str(self.stop_index)+")"
def read_next(last_line_read):
    file_object = open(sys.argv[1],"r")
    seq = Sequence()
    seq.last_line = last_line_read
    READING = False
    for i in range(last_line_read):
        file_object.readline()
    while True:
        line = file_object.readline()
        if not line:
            seq.last_line = -1
            break
        if(line[0]==">" and not READING): # if greater than symbol, line is header, starting new sequence
            #strip > sign
            line = line[1:]
            #strip newline character
            line = line[:len(line)-1]
            #split by space
            split_line = line.split(" ")
            #get name
            seq.id = split_line[0]
            #get mutation indices
            seq.mutation_indices = split_line[1:]
            if(len(seq.mutation_indices)!=0):
                seq.mutation_indices = [int(i) for i in seq.mutation_indices[0].split(",")]
            #remove \n from last mutation index
            READING = True
        elif(line[0]==">" and READING):#next sequence, stop
            break
        elif(READING):
            seq.body=seq.body+line
        seq.last_line = seq.last_line +1
    file_object.close
    seq.body = seq.body.replace("\n","")
    seq.length = len(seq.body)
    return seq
def compare_sequences(parent, seq):
    mutation_indices = []
    length = len(parent) #without loss of generality
    if(length != len(seq)):#sequence length not equal
        return
    if(parent==seq):#if the sequences are equal, we are done
        return mutation_indices
    for i in range(length):#otherwise, loop and comapre each character
        if(parent[i]!=seq[i]):
            mutation_indices.append(i)
    return mutation_indices
def recursive_read_fasta(last_line=0):
    if(last_line)==-1:
        return
    seq = read_next(last_line)
    recursive_read_fasta(seq.last_line)
def scan_for_query(seq):
    query_matches = []
    for i in range(seq.length-len(searchpattern)):
        if(seq.body[i:i+len(searchpattern)]==searchpattern):
            query_matches.append(i)
    return query_matches
def scan_for_gene(seq):
    START_CODONS = ["ATG"]
    STOP_CODONS = ["TAA","TAG","TGA"]
    WAITING_FOR_START = -1
    WAITING_FOR_STOP = -2
    state = WAITING_FOR_START
    genes = []
    gene = Gene()
    # Open Reading Frame at index 0
    for i in range(0,seq.length,3):
        current_codon = seq.body[i:i+3]
        if(current_codon in START_CODONS or current_codon in STOP_CODONS):
            print(current_codon,i)
        # if(state == WAITING_FOR_START):
        #     if(current_codon in START_CODONS): #found
        #         print(current_codon)
        #         print(i)
        # elif(state == WAITING_FOR_STOP):
        #     if(current_codon in STOP_CODONS): #found
        #         gene.stop_index=i
    return genes
def write_seq_to_file(seq,file_object):
    file_object.write(str(seq)+"\n")
def main():
    start_time = time.time()
    # parent = read_next(0)
    # last_line = parent.last_line
    # mutation_indices = compare_sequences(parent.body,parent.body)
    # if(parent.mutation_indices==mutation_indices):
    #     parent.recovery = True
    # else:
    #     parent.recovery = False
    # parent.query_matches = (scan_for_query(parent))
    # file_object = open(outfile,"w+")
    # write_seq_to_file(parent,file_object)
    # file_object.close
    # file_object = open(outfile,"a")
    # while(last_line!=-1):
    #     seq = read_next(last_line)
    #     last_line = seq.last_line
    #     mutation_indices = compare_sequences(parent.body,seq.body)
    #     if(seq.mutation_indices==mutation_indices):
    #         seq.recovery = True
    #     else:
    #         seq.recovery = False
    #     seq.query_matches = scan_for_query(seq)
    #     write_seq_to_file(seq,file_object)
    # file_object.close
    seq = read_next(0)
    genes = scan_for_gene(seq)
    for gene in genes:
        print(str(gene))
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    main()
