# author: Andrew Smith
# date: 100420 @ 17:43
# file: structures.py
# description: contains sequence and gene objects
class Sequence:
    def __init__(self,offset=0):
        self.id = "NONE"
        self.mutation_indices = 0
        self.length = 0
        self.body = ""
        self.offset = offset
        self.length = 0
        self.recovery = True
        self.query_matches = 0
        self.genes = 0
    def __str__(self):
        ret = self.id+", "+str(self.length)+", "+str(self.recovery)+", "+str(self.query_matches)
        for gene in self.genes:
            ret = ret + ", "+str(gene)
        return ret
class Gene:
    def __init__(self,start_index=0):
        self.start_index = start_index
        self.stop_index = 0
    def __str__(self):
        return "("+str(self.start_index)+","+str(self.stop_index)+")"
    def __repr__(self):
        return str(self)