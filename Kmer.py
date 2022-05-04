class Kmer:

    def __init__(self, seq):
        self.seq = seq
        self.used = False
        self.times_seen = 1

    def increment_seen(self):
        self.times_seen += 1

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.seq == other.seq
        else:
            return False
    
        