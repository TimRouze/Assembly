#include "assembler.h"

<BloomFilter, kmer> struct pair construct_filter(string filename, int k, int s){
    ifstream fasta(filename);
    BloomFilter * b_filter = new BloomFilter(k, s);
    string sequence = "";
	if(fasta.good()){
		string desc;
		while(not fasta.eof()){
			getline(fasta, desc);
			getline(fasta, line);
            for(int i = 0; i < line.length() - k + 1; i++){
                kmer k_mer = str2num(substr(line, i, k));
                if(i == 0){
                    kmer first = k_mer;
                }
                if (!b_filter->contains(k_mer)){
                    b_filter->insert(k_mer);
                }
            }
        }
    }
    return pair(b_filter, first);
}

string create_contig(string filename, int k, int h, int s){
    pair<BloomFilter, kmer> res = read_file(filename, k, s);
    kmer first = res.second;
    BloomFilter index = res.first;

    vector<bool> succs = successor(curr_kmer, index);
    if(accumulate(succs.begin(), succs.end(), 0) == 1){
        for(int i = 0; i <= 3; i++){
            if(succs[i] == 1){
                contig += "ACTG"[i];
            }
        }
    }
    return contig
}