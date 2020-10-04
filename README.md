# fasta utilities

# main.py
Outcomes:

Develop functional familiarity with Python and its applications in Bioinformatics
Develop introductory concepts in string pattern matching
Develop a functional understanding of computational and space complexity of algorithms in Bioinformatics
Overall description:

Write a Python program that will analyze the contents of a large fasta file and perform string matching. 

Details:

Your program will require 3 command-line arguments as described below: 

infile: The name of the input sequence file in fasta format.
outfile: The output file to which the results will be written. 
searchpattern: A valid genomic inquiry sequence. 

Each sequence in the input fasta file will contain an ID, a sequence of command delimited numbers indicating the points of mutation, followed by the sequence of the gene. 
The first sequence in the file will be used as a template to which all other sequences will be compared. 
Your program will produce the following information for each sequence in the fasta file (including the parent sequence):
The identity of the gene
The gene length
Ture/False indicating the successful recovery of the points of mutation
The start and stop point of any present gene in this sequence. The positions (start or stop) should be reported by assuming the first character of the sequence is at location 0. Each pair start/stop pair should be reported as a tuple. For instance, the start position of 12 and the endpoint of 892 should be reported as (12, 892).
The starting position of the search query (if present) should be reported in square brackets. For example, the appearance of the query sequence at position 1092 should be reported as [1092].
All output elements should be separated by commas. 
Your program should report your program execution time before its termination.
Example:

Executed command: “analyzeSequences.py infile outfile "TCCG"” with the following input file will produce the output file below:

Input file :

>S0
AATTCCGGAA
>S1 5
AATTCAGGAA
>S2 2
AACTCCGGAA
Output file:
S0, 10, True, [3]
S1, 10, True
S2, 10, True, [3]
Note that the above example contains no genes. 

# utilities/generateSequences.py
Outcomes:

Develop functional familiarity with Python and its applications in Bioinformatics
Develop knowledge of biological sequences
Develop a familiarity with format and mechanism of storing biological sequences
Overall description:

Write a Python program that will generate a random gene sequence of size l bases specified by the user. In addition, the program will generate m derivatives of the parent sequence by mutating p% of the bases in the parent gene. The final m+1 derivative genes will then be stored in a fasta file.

Details:

Your program will receive 4 items from the user in the form of a command line argument. The four parameters are as follows:

l: The length of the parent sequence.
m: The number of derivative sequences to be generated from the parent sequence. 
p: the percentage of mutations each child sequence should have when compared to the parent (rounded
down).
Output file name : name of the fasta output file. 

Using these parameters, your program will create m+1 sequences of length l with p% mutations compared to the parent sequence recorded in the fasta formatted file name specified by the user. The parent sequence will always be the first entry in the output file and it will be named “S0.” All subsequent children sequences will be named as Si, where i is between 1 and m. Following the sequence name, the location of individual mutations for each sequence should be listed in a comma-separated format.

Example:

Executed command: “generateSequences.py 10 3 1 test.fasta”

Output file named test.fasta with the following content:

>S0
AATTCCGGAA
>S1 5
AATTCAGGAA
>S2 2
AACTCCGGAA