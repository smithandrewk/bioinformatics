# author: Andrew Smith
# date: 100420 @ 17:45
# file: fasta_utilities.py
# description: to perform useful functions in relation to analyzing RNA from fasta
import sys
sys.path.append("..")
from utils.structures import *

def get_next_sequence(infile,outfile,file_mode,searchpattern,offset=0,parent=0):
    """
    Called for each sequence in fasta file, where interesting attributes are derived such as:
    - if the mutations were successfully recovered from the fasta file
    - the locations of query matches for given query in command line
    - start and end points of genes present in RNA sequence
    Finally, writes interesting analytics to output file
    """
    seq = read_next_sequence(infile,offset)
    if(parent==0):
        seq.recovery = recovered_mutations(seq,seq)
    else:
        seq.recovery = recovered_mutations(parent,seq)
    seq.query_matches = scan_for_query(seq,searchpattern)
    seq.genes = scan_for_genes(seq)
    seq.genes.sort(key=lambda x: x.start_index)
    write_seq_to_file(seq,outfile,file_mode)
    return seq
def read_next_sequence(infile,offset):
    """
    Called for each sequence to parse fasta and to read header and body of sequence into memory
    """
    file_object = open(infile,"r")
    seq = Sequence(offset)
    READING_BODY = False
    #skip to next sequence
    file_object.seek(offset)
    while True:
        line = file_object.readline()
        if not line:
            seq.offset = -1
            break
        if(line[0]==">" and not READING_BODY): # if greater than symbol, line is header, starting new sequence
            line = line[1:]#strip > sign
            line = line[:len(line)-1]#strip newline character
            split_line = line.split(" ")#split by space
            seq.id = split_line[0]#get name
            seq.mutation_indices = split_line[1:]#get mutation indices
            if(len(seq.mutation_indices)!=0): #correction for parsing string
                seq.mutation_indices = [int(i) for i in seq.mutation_indices[0].split(",")]
            READING_BODY = True 
        elif(line[0]==">" and READING_BODY):#next sequence, stop
            break #sequence complete
        elif(READING_BODY):
            seq.body=seq.body+line
        seq.offset = seq.offset +len(line)
    file_object.close
    seq.body = seq.body.replace("\n","")
    seq.length = len(seq.body)
    return seq
def recovered_mutations(parent, seq):
    """
    To analyze whether mutations given in fasta file are the same as mutations in sequence in memory
    """
    mutation_indices = []
    length = len(parent.body) #without loss of generality
    if(length != len(seq.body)):#sequence length not equal
        return
    for i in range(length):#otherwise, loop and comapre each character
        if(parent.body[i]!=seq.body[i]):
            mutation_indices.append(i)
    if(mutation_indices==seq.mutation_indices):
        return True
    else:
        return False
def scan_for_query(seq,searchpattern):
    """
    Simple windowed string matching for a string and a given search pattern.
    Returns list of locations of query matches.
    """
    query_matches = []
    for i in range(seq.length-len(searchpattern)):
        if(seq.body[i:i+len(searchpattern)]==searchpattern):
            query_matches.append(i)
    return query_matches
def write_seq_to_file(seq,outfile,mode):
    """
    Wrapper to write sequence to file object with a certain specified file mode
    """
    file_object = open(outfile,mode)
    file_object.write(str(seq)+"\n")
    file_object.close
def scan_for_genes(seq):
    """
    To call get_genes_from_rf for each reading frame, then return the list of genes
    """
    genes = []
    for rf in range(3):
        genes = genes+get_genes_from_rf(seq,rf)
    return genes
def get_genes_from_rf(seq,rf):
    """
    To obtain start and stop location of all genes within a reading frame
    """
    START_CODONS = ["AUG"]
    STOP_CODONS = ["UAA","UAG","UGA"]
    curr_genes  = []
    genes = []
    # Open Reading Frame at index 0
    for i in range(rf,seq.length,3):
        current_codon = seq.body[i:i+3]
        if(current_codon in START_CODONS):
            curr_genes.append(Gene(start_index=i))
        if(current_codon in STOP_CODONS):
            for gene in curr_genes:
                gene.stop_index=i+3
                genes.append(gene)
            curr_genes = []
    return genes
