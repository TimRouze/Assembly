import sys, argparse
from Bio import SeqIO
from Kmer import *
from mimetypes import guess_type
from functools import partial
import time, psutil, doctest, gzip

def construct_index(filename, k):
    kmer_index = {}
    encoding = guess_type(filename)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(filename) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            curr_seq = record.seq
            for i in range(len(curr_seq) - k+1):
                if curr_seq[i:i+k] in kmer_index:
                    kmer_index[curr_seq[i:i+k]].increment_seen()
                else:
                    kmer_index[curr_seq[i:i+k]] = Kmer(curr_seq[i:i+k])
    to_del = []
    for kmer in kmer_index.keys():
        if kmer_index[kmer].times_seen <= 3:
            to_del.append(kmer)
    for elem in to_del:
        del kmer_index[elem]
    return kmer_index

def successor(kmer, index):
    succs = [0, 0, 0, 0]
    next = 'ACTG'
    kmer = kmer[1:]
    for i in range(0,4):
        if kmer+next[i] in index:
            succs[i] = 1
    return succs

def predecessor(kmer, index):
    preds = [0, 0, 0, 0]
    next = 'ACTG'
    kmer = kmer[0:len(kmer)-1]
    for i in range(0,4):
        if next[i]+kmer in index:
            preds[i] = 1
    return preds

def build_forwards(curr_kmer, index, k):
    cpt = 0
    curr_contig = curr_kmer
    while True:
        succs = successor(curr_kmer, index)
        if sum(succs) == 1:
            for i in range(0,4):
                if succs[i] == 1:
                    if not index[curr_contig[cpt+1:]+"ACTG"[i]].used:
                        index[curr_contig[cpt+1:]+"ACTG"[i]].used = True
                        curr_contig += "ACTG"[i]
                        cpt += 1
                        break
                    else:
                        print("forward kmer already used: " + curr_kmer)
                        return curr_contig
            curr_kmer = curr_contig[cpt:]
        elif sum(succs) >= 1:
            print("several paths")
            return curr_contig
        else:
            print("dead end")
            return curr_contig

def build_backwards(curr_kmer, index, k):
    start_kmer = curr_kmer
    curr_contig = curr_kmer
    while True:
        preds = predecessor(curr_kmer, index)
        if sum(preds) == 1:
            for i in range(0,4):
                if preds[i] == 1: 
                    if not index["ACTG"[i] + curr_contig[0:k-1]].used:
                        index["ACTG"[i] + curr_contig[0:k-1]].used = True
                        curr_contig = "ACTG"[i] + curr_contig
                        break
                    else:
                        print("backward kmer already used: " + curr_kmer)
                        return curr_contig
            curr_kmer = curr_contig[0:k]
        elif sum(preds) >= 1:
            print("several paths")
            return curr_contig
        else:
            print("dead end")
            return curr_contig

def build_contig(curr_kmer, index, k):
    back_contig = build_backwards(curr_kmer, index, k)
    forward_contig = build_forwards(curr_kmer, index, k)
    return back_contig + forward_contig[k:]

def assemble(filename, k, output_file):
    res = []
    index = construct_index(filename, k)
    n = 1
    with open(output_file, 'w') as f:
        for elem in index.keys():
            if not index[elem].used:
                print(f"> contig n°{n}", file = f)
                n += 1
                print(build_contig(elem, index, k), file = f)

if __name__ == "__main__":
    #Command line options
    parser = argparse.ArgumentParser(description = 'DNA assembly tool.')
    parser.add_argument('reads', help='Fasta file containing reference sequence')
    parser.add_argument('out', help='Output file')
    parser.add_argument('k', type=int, help='Word size of the kmers the assembler should use')
    args = parser.parse_args(sys.argv[1:])
    assemble(args.reads, args.k, args.out)