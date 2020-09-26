#!/usr/bin/env python3
#author: Andrew Smith
#last edited 092420 @ 10:37
#description: to read in a fasta file

import sys
import time

#Proper Script Usage
if(len(sys.argv)!=2):
    print("Usage: ./read_fasta.py {filename}")
    sys.exit(0)
class Sequence:
    def __init__(self):
        self.id = "NONE"
        self.mutation_indices = 0
        self.length = 0
        self.body = ""
        self.last_line = 0
        self.length = 0
    def __str__(self):
        if(len(self.mutation_indices)==0):
            return self.id+self.body
        return self.id+" "+"".join(self.mutation_indices)+self.body
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
            #split by space
            split_line = line.split(" ")
            #get name
            seq.id = split_line[0]
            #get mutation indices
            seq.mutation_indices = split_line[1:]
            #remove \n from last mutation index
            READING = True
        elif(line[0]==">" and READING):#next sequence, stop
            break
        elif(READING):
            seq.body=seq.body+line
        seq.last_line = seq.last_line +1
    file_object.close
    seq.length = len(seq.body)
    return seq
def read_fasta(last_line=0):
    if(last_line)==-1:
        return
    seq = read_next(last_line)
    print(seq,end="")
    read_fasta(seq.last_line)
def main():
    start_time = time.time()
    read_fasta()
    print("--- %s seconds ---" % (time.time() - start_time))
if __name__ == '__main__':
    main()
