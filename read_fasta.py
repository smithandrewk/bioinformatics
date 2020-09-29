from structures import *
def get_next_sequence(infile,outfile,file_mode,searchpattern,last_line_read=0,parent=0):
    seq = read_next_sequence(infile,last_line_read)
    if(parent==0):
        seq.recovery = recovered_mutations(seq,seq)
    else:
        seq.recovery = recovered_mutations(parent,seq)
    seq.query_matches = scan_for_query(seq,searchpattern)
    write_seq_to_file(seq,outfile,file_mode)
    return seq
def read_next_sequence(infile,last_line_read):
    file_object = open(infile,"r")
    seq = Sequence(last_line_read)
    READING_BODY = False
    #skip to next sequence
    file_object.seek(last_line_read)
    while True:
        line = file_object.readline()
        if not line:
            seq.last_line = -1
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
        seq.last_line = seq.last_line +len(line)
    file_object.close
    seq.body = seq.body.replace("\n","")
    seq.length = len(seq.body)
    return seq
def recovered_mutations(parent, seq):
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
    query_matches = []
    for i in range(seq.length-len(searchpattern)):
        if(seq.body[i:i+len(searchpattern)]==searchpattern):
            query_matches.append(i)
    return query_matches
def write_seq_to_file(seq,outfile,mode):
    file_object = open(outfile,mode)
    file_object.write(str(seq)+"\n")
    file_object.close
def scan_for_genes(seq):
    START_CODONS = ["ATG"]
    STOP_CODONS = ["TAA","TAG","TGA"]
    WAITING_FOR_START = -1
    WAITING_FOR_STOP = -2
    state = WAITING_FOR_START
    genes = []
    # Open Reading Frame at index 0
    for i in range(0,seq.length,3):
        current_codon = seq.body[i:i+3]
        if(current_codon in START_CODONS or current_codon in STOP_CODONS):
            print(current_codon,i)
    return genes