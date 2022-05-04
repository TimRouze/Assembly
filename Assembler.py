import sys, argparse
from Bio import SeqIO
from Kmer import *
from mimetypes import guess_type
from functools import partial
import gzip

TIP = 0
FORK = 1
SIMPLE_PATH = 2
BUBBLE = 3

def construct_index(filename, k):
    """
    Build the index of the K-mers present in a set of reads.

    Parameters
    ----------
    filename : str
        File description of a fasta file containing a set of reads.
    k : int
        The word size of the K-mers

    Returns
    -------
    dict of {str : Kmer.Kmer}
        A dictionary with key the sequence of a K-mer and with value the K-mer object itself.
    """
    kmer_index = {}
    # We open the .fasta or .gz file and use the functions provided by biopython to read the file
    encoding = guess_type(filename)[1]
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
    with _open(filename) as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            curr_seq = record.seq
            # We fetch each K-mer of word size k
            for i in range(len(curr_seq) - k+1):
                # If the K-mer is already present in the index, we increment it's occurence
                if curr_seq[i:i+k] in kmer_index:
                    kmer_index[curr_seq[i:i+k]].increment_seen()
                # If the K-mer isn't present in the index, we insert it
                else:
                    kmer_index[curr_seq[i:i+k]] = Kmer(curr_seq[i:i+k])

    # We filter out the K-mers present less than 3 times in the reads
    to_del = []
    for kmer in kmer_index.keys():
        if kmer_index[kmer].times_seen <= 3:
            to_del.append(kmer)
    for elem in to_del:
        del kmer_index[elem]
    return kmer_index

def successor(kmer, index):
    """
    Determine the successors of a given K-mer present in the index.

    Parameters
    ----------
    kmer : str
        The sequence of the Kmer whose succesors we are looking for.
    index : dict of {str : Kmer.Kmer}
        The index to look for succesors in.

    Returns
    -------
    list of str
        A list of all the successors of kmer present in the index, and which have not already been used to build a contig.
    """
    succs = []
    next = 'ACTG'
    kmer = kmer[1:]
    for i in range(0,4):
        if kmer+next[i] in index and not index[kmer+next[i]].used:
            succs.append(kmer+next[i])
    return succs

def predecessor(kmer, index):
    """
     Determine the predecessors of a given K-mer present in the index.

    Parameters
    ----------
    kmer : str
        The sequence of the Kmer whose predecessors we are looking for.
    index : dict of {str : Kmer.Kmer}
        The index to look for predecessors in.

    Returns
    -------
    list of str
        A list of all the predecessors of kmer present in the index, and which have not already been used to build a contig.
    """
    preds = []
    next = 'ACTG'
    kmer = kmer[0:len(kmer)-1]
    for i in range(0,4):
        if next[i]+kmer in index and not index[next[i]+kmer].used:
            preds.append(next[i]+kmer)
    return preds

def build_forwards(curr_contig, curr_kmer, index, k):
    """
    Build a contig using a start K-mer and the index in the direction left to right.

    Parameters
    ----------
    curr_contig : str
        The contig built by build_backwards from the same kmer.
    curr_kmer : str
        The sequence of the K-mer to start building the contig from.
    index : dict of {str : Kmer.Kmer}
        The index of K-mers.
    k : int
        The word size of K-mers.

    Returns
    -------
    list of str
        A list containing the contigs built.
    """
    # We continue building while we can
    while True:
        succs = successor(curr_kmer, index)
        # If there is no fork in the De Brujin graph, we simply construc a unitig bp by bp
        if len(succs) == 1:
            index[succs[0]].used = True
            curr_contig += succs[0][-1]
            curr_kmer = succs[0]
        # If there is a simple fork, we need to check whether it is a pattern we can try and simplify
        elif len(succs) == 2:
            print("two paths forwards")
            contigs = []
            codes = []
            # We look at each side of the fork to see if we have a tip/bubble pattern
            for elem in succs:
                contig, code = second_extension_forward(index, elem, k)
                contigs.append(contig)
                codes.append(code)
            # We use the check_codes_forward function to delete and/or merge some contigs if we have a tip/bubble pattern
            return check_codes_forward(contigs, codes, k, curr_contig, index)
        # If we get to a dead end in the De Brujin graph or a fork with 3 or 4 successors, we cannot build anymore.
        else:
            print("dead end or too many paths forwards")
            return [curr_contig]

def build_backwards(curr_kmer, index, k):
    """
    Build a contig using a start K-mer and the index in the direction right to left.

    Parameters
    ----------
    curr_kmer : str
        The sequence of the K-mer to start building the contig from.
    index : dict of {str : Kmer.Kmer}
        The index of K-mers.
    k : int
        The word size of K-mers.

    Returns
    -------
    list of str
        A list containing the contigs built.
    """
    curr_contig = curr_kmer
    # We continue building while we can
    while True:
        preds = predecessor(curr_kmer, index)
        # If there is no fork in the De Brujin graph, we simply construc a unitig bp by bp
        if len(preds) == 1:
            index[preds[0]].used = True
            curr_contig = preds[0][0] + curr_contig
            curr_kmer = preds[0]
        # If there is a simple fork, we need to check whether it is a pattern we can try and simplify
        elif len(preds) == 2:
            print("two paths backwards")
            contigs = []
            codes = []
            # We look at each side of the fork to see if we have a tip/bubble pattern
            for elem in preds:
                contig, code = second_extension_backward(index, elem, k)
                contigs.append(contig)
                codes.append(code)
            # We use the check_codes_backward function to delete and/or merge some contigs if we have a tip/bubble pattern
            return check_codes_backward(contigs, codes, k, curr_contig, index)
        # If we get to a dead end in the De Brujin graph or a fork with 3 or 4 predecessors, we cannot build anymore.
        else:
            print("dead end or too many paths backwards")
            return [curr_contig]

def second_extension_backward(index, curr_kmer, k):
    """
    Build a contig after a simple fork in the direction right to left.

    Parameters
    ----------
    index : dict of {str : Kmer.Kmer}
        The index of K-mers.
    curr_kmer : str
        The sequence of the K-mer to start building the contig from.
    k : int
        The word size of K-mers.

    Returns
    -------
    tuple of (str, int)
        A pair containing the sequence of the contig built, and a code corresponding to the reason why the extension stopped.
    """
    curr_contig = curr_kmer
    # We continue building while we can
    while True:
        preds = predecessor(curr_kmer, index)
        # If there is no fork in the De Brujin graph, we simply construct a unitig bp by bp
        if len(preds) == 1:
            index[preds[0]].used = True
            curr_contig = preds[0][0] + curr_contig
            curr_kmer = preds[0]
        # If there is a fork, we cannot continue extending, and we return the contig and the FORK exit code
        elif len(preds) >= 1:
            print("more than one path, pattern is too complex.")
            return (curr_contig, FORK)
        # If the contig comes to a dead end, and is shortand is less than 3 times k long, we return the contig and the TIP exit code
        elif len(preds) == 0 and len(curr_contig) <= 3*k:
            print("This is a tip")
            return (curr_contig, TIP)
        # If the contig comes to a dead end, but is not short enough to be a tip, we return the contig and the SIMPLE_PATH exit code
        else:
            print("Simple path (plan lol)")
            return (curr_contig, SIMPLE_PATH)

def second_extension_forward(index, curr_kmer, k):
    """
    Build a contig after a simple fork in the direction left to right.

    Parameters
    ----------
    index : dict of {str : Kmer.Kmer}
        The index of K-mers.
    curr_kmer : str
        The sequence of the K-mer to start building the contig from.
    k : int
        The word size of K-mers.

    Returns
    -------
    tuple of (str, int)
        A pair containing the sequence of the contig built, and a code corresponding to the reason why the extension stopped.
    """
    curr_contig = curr_kmer
    # We continue building while we can
    while True:
        succs = successor(curr_kmer, index)
        # If there is no fork in the De Brujin graph, we simply construct a unitig bp by bp
        if len(succs) == 1:
            index[succs[0]].used = True
            curr_contig += succs[0][-1]
            curr_kmer = succs[0]
        # If there is a fork, we cannot continue extending, and we return the contig and the FORK exit code
        elif len(succs) >= 1:
            print("more than one path, pattern is too complex.")
            return (curr_contig, FORK)
        # If the contig comes to a dead end, and is less than 3 times k long, we return the contig and the TIP exit code
        elif len(succs) == 0 and len(curr_contig) <= 3*k:
            print("This is a tip")
            return (curr_contig, TIP)
        # If the contig comes to a dead end, but is not short enough to be a tip, we return the contig and the SIMPLE_PATH exit code
        else:
            print("Simple path (plan lol)")
            return (curr_contig, SIMPLE_PATH)

def check_codes_forward(contigs, codes, k, curr_contig, index):
    """
    Delete and/or merge some contigs if they correspond to tip/bubble patterns in the direction left to right.

    Parameters
    ----------
    contigs : list of str
        List containing two contigs potentially corresponding to a tip/bubble pattern.
    codes : list of int
        List of exit codes provided by second_extension_forward or second_extension_backward.
    curr_contig : str
        The contig from which the tip/bubble pattern extends from.
    index : dict of {str : Kmer.Kmer}
        The index of K-mers.
    k : int
        The word size of K-mers.

    Returns
    -------
    list of str
        Initial list of contigs after deleting and/or merging some contigs depending on the codes provided.
    """
    # If one side of the fork contains a short dead end, and the other doesn't, this is a tip pattern, we delete the short dead end
    # and merge the remaining side of the fork with the source contig
    if codes[0] ==  TIP and codes[1] != TIP:
        curr_contig += contigs[1][k-1:]
        return [curr_contig]
    elif codes[0] != TIP and codes[1] == TIP:
        curr_contig += contigs[0][k-1:]
        return [curr_contig]
    # If both sides of the fork contain a simple path, are of same length and have the same ending K-mer, this is a bubble pattern, we 
    # must check the scores of each side of the K-mer
    elif codes[0] == SIMPLE_PATH and codes[1] == SIMPLE_PATH and len(contigs[0]) == len(contigs[1]) and contigs[0][-k+1] == contigs[0][-k+1]:
        score1 = score_contig(contigs[0], k, index)
        score2 = score_contig(contigs[1], k, index)
        # If both sides have the same score, this bubble is probably due to two different haplotypes in the sequencing reads, we cannot 
        # delete one of the sides
        # If there is a disparity between the scores, the bubble is might be due to a sequencing error in the reads, we thus set a 1/5
        # threshold ratio to be sure that we are deleting a sequencing error
        if score1 == score2:
            return[curr_contig] + contigs
        else:
            s_min =  min(score1, score2)
            s_max = max(score1, score2)
            if s_min < s_max*0.2:
                if s_max == score2:
                    return [curr_contig + contigs[1][k-1:]]
                else:
                    return [curr_contig + contigs[0][k-1:]]
            else:
                return[curr_contig] + contigs
    else:
        return [curr_contig] + contigs

def check_codes_backward(contigs, codes, k, curr_contig, index):
    """
    Delete and/or merge some contigs if they correspond to tip/bubble patterns in the direction left to right.

    Parameters
    ----------
    contigs : list of str
        List containing two contigs potentially corresponding to a tip/bubble pattern.
    codes : list of int
        List of exit codes provided by second_extension_forward or second_extension_backward.
    curr_contig : str
        The contig from which the tip/bubble pattern extends from.
    index : dict of {str : Kmer.Kmer}
        The index of K-mers.
    k : int
        The word size of K-mers.

    Returns
    -------
    list of str
        Initial list of contigs after deleting and/or merging some contigs depending on the codes provided.
    """
    # If one side of the fork contains a short dead end, and the other doesn't, this is a tip pattern, we delete the short dead end
    # and merge the remaining side of the fork with the source contig
    if codes[0] ==  TIP and codes[1] != TIP:
        curr_contig = contigs[1][:-k+1] + curr_contig
        return [curr_contig]
    elif codes[0] != TIP and codes[1] == TIP:
        curr_contig = contigs[0][:-k+1] + curr_contig
        return [curr_contig]
    # If both sides of the fork contain a simple path, are of same length and have the same ending K-mer, this is a bubble pattern, we 
    # must check the scores of each side of the K-mer
    elif codes[0] == SIMPLE_PATH and codes[1] == SIMPLE_PATH and len(contigs[0]) == len(contigs[1]) and contigs[0][-k+1] == contigs[0][-k+1]:
        score1 = score_contig(contigs[0], k, index)
        score2 = score_contig(contigs[1], k, index)
        # If both sides have the same score, this bubble is probably due to two different haplotypes in the sequencing reads, we cannot 
        # delete one of the sides
        # If there is a disparity between the scores, the bubble is might be due to a sequencing error in the reads, we thus set a 1/5
        # threshold ratio to be sure that we are deleting a sequencing error
        if score1 == score2:
            return[curr_contig] + contigs
        else:
            s_min =  min(score1, score2)
            s_max = max(score1, score2)
            if s_min < s_max*0.2:
                if s_max == score2:
                    return [contigs[1][:-k+1] + curr_contig]
                else:
                    return [contigs[0][:-k+1] + curr_contig]
            else:
                return[curr_contig] + contigs
    else:
        return [curr_contig] + contigs

def score_contig(contig, k, index):
    """
    Assign a score to a contig corresponding to the sum of the occurences of the Kmers composing the contig.

    Parameters
    ----------
    contig : str
        sequence of the contig we wish to score.
    k : int
        The word size of Kmers.
    index : dict of {str : Kmer.Kmer}
        The index of K-mers.

    Returns
    -------
    int
        The score of the contig.
    """
    score = 0
    for i in range(0, len(contig)-k):
        score += index[contig[i:i+k]].times_seen
    return score

def build_contigs(index, k):
    """
    Generate all the contigs using a De Brujin graph. 

    Parameters
    ----------
    index : dict of {str : Kmer.Kmer}
        The index of K-mers.
    k : int
        The word size of K-mers.

    Returns
    -------
    list of str
        List of all the contigs generated.
    """
    results = []
    for kmer in index.keys():
        if not index[kmer].used:
            back_contig = build_backwards(kmer, index, k)
            forward_contig = build_forwards(back_contig[0], kmer, index, k)
            if len(back_contig) > 1:
                results = results + back_contig[1:]
            if len(forward_contig) > 1:
                results = results + forward_contig[1:]
            results.append(forward_contig[0])
    return results

def assemble(filename, k, output_file):
    """
    Assemble a genome using from a set of sequencing reads.

    Parameters
    ----------
    filename : str
        File description of a fasta file containing a set of reads.
    k : int
        The word size of K-mers.
    output_file : str
        File description of an output file to write the results in.
    l : int
        Minimum number of times a K-mer mist be present in the reads for it to be kept in the index.
    """
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