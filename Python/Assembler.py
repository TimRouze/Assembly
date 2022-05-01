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
        #print(preds)
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
        if contig.seq not in contigs:
            contigs[contig.seq] = contig
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
                contigs[contig.seq].add_pred(backward_contig)
                if backward_contig.seq not in contigs:
                    contigs[backward_contig.seq] = backward_contig

def multi_build_forwards(curr_kmer, index, contigs):
    checked_kmer = [curr_kmer]
    k = len(curr_kmer)
    to_test = {}
    contig, curr_kmer = build_forwards(curr_kmer, index, 0)
    if len(contig.seq) > k:
        if contig.seq not in contigs:
            contigs[contig.seq] = contig
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
                contigs[contig.seq].add_succ(forward_contig)
                if forward_contig.seq not in contigs:
                    contigs[forward_contig.seq] = forward_contig

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
    contigs = {}
    curr_kmer = list(index.keys())[0]
    multi_build_backwards(curr_kmer, index, contigs)
    multi_build_forwards(curr_kmer, index, contigs)
    for out_key in contigs.keys():
        for in_key in contigs.keys():
            if contigs[out_key] in contigs[in_key].preds:
                contigs[out_key].add_succ(contigs[in_key])
            if contigs[out_key] in contigs[in_key].succs:
                contigs[out_key].add_pred(contigs[in_key])
            if contigs[out_key].kmer == contigs[in_key].kmer and contigs[out_key] != contigs[in_key]:
                if contigs[out_key].is_forward:
                    contigs[out_key].add_pred(contigs[in_key])
                    contigs[in_key].add_succ(contigs[out_key])
                else:
                    contigs[out_key].add_succ(contigs[in_key])
                    contigs[in_key].add_pred(contigs[out_key])
            
    return contigs

def assemble(filename, k):
    index = construct_index(filename, k)
    contigs = create_contigs(index)
    """for key in contigs.keys():
        print("seq: " + contigs[key].seq)
        print("nb successors: " + str(len(contigs[key].succs)))
        print("nb predecessors: " + str(len(contigs[key].preds)))
    """
    for key in contigs.keys():
        
        print("------------------PREDECESSORS-------------")
        for i in range(len(contigs[key].preds)):
            print(contigs[key].preds[i].seq)
        print("------------------CONTIG-------------------")
        print(contigs[key].seq)
        print("------------------SUCCESSORS-----------------")
        for i in range(len(contigs[key].succs)):
            print(contigs[key].succs[i].seq[k:])
    


if __name__ == "__main__":
    #Command line options
    parser = argparse.ArgumentParser(description = 'Seed and extend alignment tool.')
    parser.add_argument('genome', help='Fasta file containing reference sequence')
    parser.add_argument('out', help='Output file')
    parser.add_argument('k', type=int, help='Word size of the kmers the alignment should use')
    args = parser.parse_args(sys.argv[1:])
    assemble(args.genome, args.k)
    '''for key in index.keys():
        print(key, index[key])'''
    