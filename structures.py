class Sequence:
    def __init__(self,last_line=0):
        self.id = "NONE"
        self.mutation_indices = 0
        self.length = 0
        self.body = ""
        self.last_line = last_line
        self.length = 0
        self.recovery = True
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