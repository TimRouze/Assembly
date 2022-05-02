import sys, argparse
from Bio import SeqIO
from Utils import *
from Contig import *
from mimetypes import guess_type
from functools import partial
import time, psutil, doctest, gzip

def construct_index(filename, k):
    index = {}
    encoding = guess_type(filename)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(filename) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            curr_seq = record.seq
            for i in range(len(curr_seq) - k+1):
                if curr_seq[i:i+k] in index:
                    index[curr_seq[i:i+k]] += 1
                else:
                    index[str(curr_seq[i:i+k])] = 1

    return index

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

def build_forwards(curr_kmer, index, cpt):
    start_kmer = curr_kmer
    curr_contig = curr_kmer
    k = len(curr_kmer)
    while True:
        succs = successor(curr_kmer, index)
        if sum(succs) == 1:
            for i in range(0,4):
                if succs[i] == 1:
                    curr_contig += "ACTG"[i]
                    cpt += 1
                    break
            curr_kmer = curr_contig[cpt:cpt + k+1]
        elif sum(succs) >= 1:
            return (Contig(curr_contig, start_kmer, True), curr_kmer)
        else:
            return (Contig(curr_contig, start_kmer, True), None)

def build_backwards(curr_kmer, index):
    start_kmer = curr_kmer
    curr_contig = curr_kmer
    k = len(curr_kmer)
    while True:
        preds = predecessor(curr_kmer, index)
        if sum(preds) == 1:
            for i in range(0,4):
                if preds[i] == 1:
                    curr_contig = "ACTG"[i] + curr_contig
                    break
            curr_kmer = curr_contig[0:k]
        elif sum(preds) >= 1:
            return (Contig(curr_contig, start_kmer, False), curr_kmer)
        else:
            return (Contig(curr_contig, start_kmer, False), None)

def multi_build_backwards(curr_kmer, index, contigs):
    checked_kmer = [curr_kmer]
    k = len(curr_kmer)
    to_test = {}
    contig, curr_kmer = build_backwards(curr_kmer, index)
    if len(contig.seq) > k:
        if contig not in contigs:
            contigs.append(contig)
        if curr_kmer != None:
            add_preds_to_check(to_test, contig, curr_kmer, index)
    for contig in to_test.keys():
        for kmer in to_test[contig]:
            if kmer not in checked_kmer:
                checked_kmer.append(kmer)
                backward_contig, curr_kmer = build_backwards(kmer, index)
                if curr_kmer != None and curr_kmer not in checked_kmer:
                    multi_build_backwards(curr_kmer, index, contigs)
                    #add_preds_to_check(to_test, backward_contig, curr_kmer, index)
                contigs[contigs.index(contig)].add_pred(backward_contig)
                if backward_contig.seq not in contigs:
                    contigs.append(backward_contig)

def multi_build_forwards(curr_kmer, index, contigs):
    checked_kmer = [curr_kmer]
    k = len(curr_kmer)
    to_test = {}
    contig, curr_kmer = build_forwards(curr_kmer, index, 0)
    if len(contig.seq) > k:
        if contig.seq not in contigs:
            contigs.append(contig)
        if curr_kmer != None:
            add_succs_to_check(to_test, contig, curr_kmer, index)
    for contig in to_test.keys():
        for kmer in to_test[contig]:
            if kmer not in checked_kmer:
                checked_kmer.append(kmer)
                forward_contig, curr_kmer = build_forwards(kmer, index, 0)
                if curr_kmer != None and curr_kmer not in checked_kmer:
                    multi_build_forwards(curr_kmer, index, contigs)
                    #add_succs_to_check(to_test, forward_contig, curr_kmer, index)
                contigs[contigs.index(contig)].add_succ(forward_contig)
                if forward_contig not in contigs:
                    contigs.append(forward_contig)

def add_succs_to_check(to_test, contig, kmer, index):
    succs = successor(kmer, index)
    for i in range(0,4):
        if succs[i] == 1:
            if contig not in to_test:
                to_test[contig] = [kmer[1:]+"ACTG"[i]]
            else:
                to_test[contig].append(kmer[1:]+"ACTG"[i])

def add_preds_to_check(to_test, contig, kmer, index):
    succs = predecessor(kmer, index)
    for i in range(0,4):
        if succs[i] == 1:
            if contig not in to_test:
                to_test[contig] = ["ACTG"[i] + kmer[:len(kmer)-1]]
            else:
                to_test[contig].append("ACTG"[i] + kmer[:len(kmer)-1])

def create_contigs(index):
    contigs = []
    curr_kmer = list(index.keys())[0]
    multi_build_backwards(curr_kmer, index, contigs)
    multi_build_forwards(curr_kmer, index, contigs)
    for first in contigs:
        for second in contigs:
            if first in second.preds:
                first.add_succ(second)
            if first in second.succs:
                first.add_pred(second)
            if first.kmer == second.kmer and first != second:
                if first.is_forward:
                    first.add_pred(second)
                    second.add_succ(first)
                else:
                    first.add_succ(second)
                    second.add_pred(first)
            
    return contigs

def assemble(filename, k):
    index = construct_index(filename, k)
    contigs = create_contigs(index)
    if len(contigs) > 2:
        contigs, cur_k = multiple_assembly(filename, contigs, k, 300)

    if len(contigs) == 2:
        LETSGO = contigs[0].seq + contigs[1].seq[cur_k:] 

    print('lol' + str(len(contigs)))

    for contig in contigs:
        
        print("------------------PREDECESSORS-------------")
        for i in range(len(contig.preds)):
            print(contig.preds[i].seq)
        print("------------------CONTIG-------------------")
        print(contig.seq)
        print("------------------SUCCESSORS-----------------")
        for i in range(len(contig.succs)):
            print(contig.succs[i].seq[k:])

    #print(LETSGO)

def multiple_assembly(filename, contigs, k, k_max):
    test = 0
    cur_k = k
    while cur_k*2 < k_max and len(contigs) > 2:
        test+=1
        print(test)
        index = construct_index(filename, cur_k)
        for contig in contigs:
            if len(contig.seq) <= cur_k:
                for i in range(len(contig.seq) - cur_k+1):
                    if contig.seq[i:i+cur_k] in index:
                        index[contig.seq[i:i+cur_k]] += 1
                    else:
                        index[str(contig.seq[i:i+cur_k])] = 1
        cur_k = cur_k*2
        contigs = create_contigs(index)
    return (contigs, cur_k)

if __name__ == "__main__":
    #Command line options
    parser = argparse.ArgumentParser(description = 'DNA assembly tool.')
    parser.add_argument('reads', help='Fasta file containing reference sequence')
    parser.add_argument('out', help='Output file')
    parser.add_argument('k', type=int, help='Word size of the kmers the assembler should use')
    args = parser.parse_args(sys.argv[1:])
    assemble(args.reads, args.k)
    '''for key in index.keys():
        print(key, index[key])'''
    