import sys, argparse
from Bio import SeqIO
from Utils import *
from Contig import *
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

def construct_index_from_contigs(contigs, k, kmer_index):
    for elem in contigs:
        curr_seq = elem.seq
        for i in range(len(curr_seq) - k+1):
            if curr_seq[i:i+k] in kmer_index:
                kmer_index[curr_seq[i:i+k]].increment_seen()
            else:
                kmer_index[curr_seq[i:i+k]] = Kmer(curr_seq[i:i+k])
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

def build_forwards(curr_kmer, index, cpt, k):
    start_kmer = curr_kmer
    curr_contig = curr_kmer
    k = len(curr_kmer)
    while True:
        succs = successor(curr_kmer, index)
        if sum(succs) == 1:
            for i in range(0,4):
                if succs[i] == 1 and not index[curr_contig[cpt:]+"ACTG"[i]].used:
                    index[curr_contig[cpt:]+"ACTG"[i]].used = True
                    curr_contig += "ACTG"[i]
                    cpt += 1
                    break
            curr_kmer = curr_contig[cpt:]
        elif sum(succs) >= 1:
            return (Contig(curr_contig, start_kmer, True), curr_kmer)
        else:
            return (Contig(curr_contig, start_kmer, True), None)

def build_backwards(curr_kmer, index, k):
    start_kmer = curr_kmer
    curr_contig = curr_kmer
    k = len(curr_kmer)
    while True:
        preds = predecessor(curr_kmer, index)
        if sum(preds) == 1:
            for i in range(0,4):
                if preds[i] == 1 and not index["ACTG"[i] + curr_contig[0:k-1]].used:
                    index["ACTG"[i] + curr_contig[0:k-1]].used = True
                    curr_contig = "ACTG"[i] + curr_contig
                    break
            curr_kmer = curr_contig[0:k]
        elif sum(preds) >= 1:
            return (Contig(curr_contig, start_kmer, False), curr_kmer)
        else:
            return (Contig(curr_contig, start_kmer, False), None)

def multi_build_backwards(curr_kmer, index, contigs, k):
    checked_kmer = [curr_kmer]
    to_test = {}
    contig, curr_kmer = build_backwards(curr_kmer, index, k)
    if len(contig.seq) > k*3:
        """if contig not in contigs:
            contigs.append(contig)"""
        if curr_kmer != None:
            add_preds_to_check(to_test, contig, curr_kmer, index)
    if len(to_test)!= 0:
        for contig in to_test.keys():
            for kmer in to_test[contig]:
                if kmer not in checked_kmer:
                    checked_kmer.append(kmer)
                    backward_contig, curr_kmer = build_backwards(kmer, index, k)
                    if curr_kmer != None and curr_kmer not in checked_kmer:
                        multi_build_backwards(curr_kmer, index, contigs)
                        #add_preds_to_check(to_test, backward_contig, curr_kmer, index)
                    contigs[contigs.index(contig)].add_pred(backward_contig)
                    if backward_contig.seq not in contigs:
                        contigs.append(backward_contig)
    else:
        return

def multi_build_forwards(curr_kmer, index, contigs, k):
    checked_kmer = [curr_kmer]
    k = len(curr_kmer)
    to_test = {}
    contig, curr_kmer = build_forwards(curr_kmer, index, 0, k)
    if len(contig.seq) > k:
        if contig.seq not in contigs:
            contigs.append(contig)
        if curr_kmer != None:
            add_succs_to_check(to_test, contig, curr_kmer, index)
    for contig in to_test.keys():
        for kmer in to_test[contig]:
            if kmer not in checked_kmer:
                checked_kmer.append(kmer)
                forward_contig, curr_kmer = build_forwards(kmer, index, 0, k)
                if curr_kmer != None and curr_kmer not in checked_kmer:
                    multi_build_forwards(curr_kmer, index, contigs, k)
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

def create_contigs(index, k):
    contigs = []
    curr_kmer = list(index.keys())[0]
    multi_build_backwards(curr_kmer, index, contigs, k)
    multi_build_forwards(curr_kmer, index, contigs, k)
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
    contigs = create_contigs(index, k)
    cur_k = k
    if len(contigs) != 2:
        contigs, cur_k = multiple_assembly(filename, contigs, k, 300)

    if len(contigs) == 2:
        print(cur_k)
        LETSGO = contigs[0].seq + contigs[1].seq[cur_k:] 
        print(LETSGO)
    elif len(contigs) == 0:
        print("pas de contigs")
    else:
        print(contigs[0].seq)
        print(len(contigs))
    


    """for contig in contigs:
        
        print("------------------PREDECESSORS-------------")
        for i in range(len(contig.preds)):
            print(contig.preds[i].seq)
        print("------------------CONTIG-------------------")
        print(contig.seq)
        print("------------------SUCCESSORS-----------------")
        for i in range(len(contig.succs)):
            print(contig.succs[i].seq[k:])"""
    
    

def multiple_assembly(filename, contigs, k, k_max):
    cur_k = k
    while cur_k*2 < k_max and len(contigs) != 2:
        cur_k = cur_k*2
        print(cur_k)
        index = construct_index(filename, cur_k)
        index = construct_index_from_contigs(contigs, cur_k, index)
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
    