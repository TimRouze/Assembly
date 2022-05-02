from Bio import SeqIO
from mimetypes import guess_type
from functools import partial
import time, types, gzip
import gzip

def parseFasta(fi):
    """Read a fasta or fasta.gz file and extract the sequence within.

    Parameters
    ----------
    fi : str
        Path of a fasta or fasta.gz file.

    Returns
    -------
    str
        The sequence in the input file.
    """
    sequence = ""
    encoding = guess_type(fi)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(fi) as fasta:
        try:
            for record in SeqIO.parse(fasta, "fasta"):
                sequence += record.seq
        except:
            print('File format incorrect or file content does not correspond to fasta.')
    return sequence

def find_kmers(seq, k):
    """Read a DNA sequence and extract all the words of size k in it.

    Parameters
    ----------
    seq : str
        A DNA sequence.
    k : int
        The word size used.

    Returns
    -------
    list of str
        A list of all the words of size k found in the sequence.
    """
    kmers = []
    for i in range(len(seq) - k+1):
        kmers.append(seq[i:i+k])
    return kmers