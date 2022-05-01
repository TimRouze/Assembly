class Contig:

    def __init__(self, seq):
        self.seq = seq
        self.preds = []
        self.succs = []
        self.kmer = ""

    def __init__(self, seq, kmer, is_forward):
        self.seq = seq
        self.kmer = kmer
        self.is_forward = is_forward
        self.preds = []
        self.succs = []
        
    
    def add_pred(self, pred):
        if type(pred) == list:
            if len(pred) != 0:
                for elem in pred:
                    if pred not in self.preds:
                        self.preds.append(pred)
        else:
            if pred not in self.preds:
                    self.preds.append(pred)
    
    def add_succ(self, succ):
        if type(succ) == list:
            if len(succ) != 0:
                for elem in succ:
                    if succ not in self.succs:
                        self.succs.append(succ)
        else:
            if succ not in self.succs:
                    self.succs.append(succ)
    
    def set_kmer(self, kmer):
        self.kmer = kmer