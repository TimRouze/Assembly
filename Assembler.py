import sys, argparse
from Bio import SeqIO
from Kmer import *
from mimetypes import guess_type
from functools import partial
import time, psutil, doctest, gzip

TIP = 0
FORK = 1
SIMPLE_PATH = 2
BUBBLE = 3

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
    succs = []
    next = 'ACTG'
    kmer = kmer[1:]
    for i in range(0,4):
        if kmer+next[i] in index and not index[kmer+next[i]].used:
            succs.append(kmer+next[i])
    return succs

def predecessor(kmer, index):
    preds = []
    next = 'ACTG'
    kmer = kmer[0:len(kmer)-1]
    for i in range(0,4):
        if next[i]+kmer in index and not index[next[i]+kmer].used:
            preds.append(next[i]+kmer)
    return preds

def build_forwards(curr_kmer, index, k):
    curr_contig = curr_kmer
    while True:
        succs = successor(curr_kmer, index)
        if len(succs) == 1:
            index[succs[0]].used = True
            curr_contig += succs[0][-1]
            curr_kmer = succs[0]
        elif len(succs) == 2:
            print("two paths forwards")
            contigs = []
            codes = []
            for elem in succs:
                contig, code = second_extension_forward(index, elem, k)
                contigs.append(contig)
                codes.append(code)
            return check_codes_forward(contigs, codes, k, curr_contig, index)
        else:
            print("dead end or too many paths forwards")
            return [curr_contig]

def build_backwards(curr_kmer, index, k):
    curr_contig = curr_kmer
    while True:
        preds = predecessor(curr_kmer, index)
        if len(preds) == 1:
            index[preds[0]].used = True
            curr_contig = preds[0][0] + curr_contig
            curr_kmer = preds[0]
        elif len(preds) == 2:
            print("two paths backwards")
            contigs = []
            codes = []
            for elem in preds:
                contig, code = second_extension_backward(index, elem, k)
                contigs.append(contig)
                codes.append(code)
            return check_codes_backward(contigs, codes, k, curr_contig, index)
        else:
            print("dead end or too many paths backwards")
            return [curr_contig]

def second_extension_backward(index, curr_kmer, k):
    curr_contig = curr_kmer
    while True:
        preds = predecessor(curr_kmer, index)
        if len(preds) == 1:
            index[preds[0]].used = True
            curr_contig = preds[0][0] + curr_contig
            curr_kmer = preds[0]
        elif len(preds) >= 1:
            print("more than one path, pattern is too complex.")
            return (curr_contig, FORK)
        elif len(preds) == 0 and len(curr_contig) <= 3*k:
            print("This is a tip")
            return (curr_contig, TIP)
        else:
            print("un autre cas (bulle? len > 3*k?)")
            return (curr_contig, SIMPLE_PATH)

def second_extension_forward(index, curr_kmer, k):
    curr_contig = curr_kmer
    while True:
        succs = successor(curr_kmer, index)
        if len(succs) == 1:
            index[succs[0]].used = True
            curr_contig += succs[0][-1]
            curr_kmer = succs[0]
        elif len(succs) >= 1:
            print("more than one path, pattern is too complex.")
            return (curr_contig, FORK)
        elif len(succs) == 0 and len(curr_contig) <= 3*k:
            print("This is a tip")
            return (curr_contig, TIP)
        else:
            print("un autre cas (bulle? len > 3*k?)")
            return (curr_contig, SIMPLE_PATH)

def check_codes_forward(contigs, codes, k, curr_contig, index):
    if codes[0] == TIP and codes[1] == TIP:
        return [curr_contig] + contigs
    elif codes[0] ==  TIP and codes[1] != TIP:
        curr_contig += contigs[1][k-1:]
        return [curr_contig]
    elif codes[0] != TIP and codes[1] == TIP:
        curr_contig += contigs[0][k-1:]
        return [curr_contig]
    elif codes[0] == SIMPLE_PATH and codes[1] == SIMPLE_PATH and len(contigs[0]) == len(contigs[1]) and contigs[0][-k+1] == contigs[0][-k+1]:
        score1 = score_contig(contigs[0], k, index)
        score2 = score_contig(contigs[1], k, index)
        
    else:
        return [curr_contig] + contigs

def check_codes_backward(contigs, codes, k, curr_contig, index):
    if codes[0] == TIP and codes[1] == TIP:
        return [curr_contig] + contigs
    elif codes[0] ==  TIP and codes[1] != TIP:
        curr_contig = contigs[1][:-k+1] + curr_contig
        return [curr_contig]
    elif codes[0] != TIP and codes[1] == TIP:
        curr_contig = contigs[0][:-k+1] + curr_contig
        return [curr_contig]
    else:
        return [curr_contig] + contigs

def score_contig(contig, k, index):
    score = 0
    for i in range(0, len(contig)-k):
        score += index[contig[i:i+k]].times_seen
    return score

def build_contigs(index, k):
    results = []
    for kmer in index.keys():
        if not index[kmer].used:
            back_contig = build_backwards(kmer, index, k)
            forward_contig = build_forwards(kmer, index, k)
            if len(forward_contig) > 1 and len(back_contig) > 1:
                results.append(back_contig[0] + forward_contig[0][k:])
                results = results + forward_contig[1:] + back_contig[1:]
            elif len(forward_contig) > 1:
                results.append(back_contig[0] + forward_contig[0][k:])
                results = results + forward_contig[1:]
            elif len(back_contig) > 1:
                results.append(back_contig[0] + forward_contig[0][k:])
                results = results + back_contig[1:]
            else:
                results.append(back_contig[0] + forward_contig[0][k:])
    return results

def assemble(filename, k, output_file):
    res = []
    index = construct_index(filename, k)
    n = 1
    with open(output_file, 'w') as f:
        res = build_contigs(index, k)
        for contig in res:
            print(f"> contig n°{n}", file = f)
            n += 1
            print(contig, file = f)

if __name__ == "__main__":
    #Command line options
    parser = argparse.ArgumentParser(description = 'DNA assembly tool.')
    parser.add_argument('reads', help='Fasta file containing reference sequence')
    parser.add_argument('out', help='Output file')
    parser.add_argument('k', type=int, help='Word size of the kmers the assembler should use')
    args = parser.parse_args(sys.argv[1:])
    assemble(args.reads, args.k, args.out)