#include "assembler.h"

pair <BloomFilter, kmer> construct_index(string filename, int k, int h, int s){
    ifstream fasta(filename);
    BloomFilter *index = new BloomFilter(h, s);
    string sequence = "";
    kmer first, k_mer;
    int n = 0;
    if (not fasta.good()) {
		cout << "Problem with file opening" << endl;
		exit(1);
	}
	else{
		string desc;
		while(not fasta.eof()){
			getline(fasta, desc);
			getline(fasta, sequence);
            for(int i = 0; i < sequence.length() - k + 1; i++){
                k_mer = str2num(sequence.substr(i, k));
                n++;
                if (!index->contains(k_mer)){
                    index->insert(k_mer);
                }
            }
        }
    }
    string file = filename.substr(filename.find_last_of("/")+1, filename.length());
    cout << file + " contains " + to_string(n) + " " + to_string(k) + "-mers."  << endl;
    cout << "The Bloom Filter caught " + to_string(index->getN()) + " unique " + to_string(k) + "-mers." << endl;

    pair<BloomFilter, kmer> res (*index, k_mer);
    return res;
}

vector<bool> successor(kmer k_mer, BloomFilter b_filter){
    vector<bool> succs(4,0);
    string next = "ACGT";
    for(int i = 0; i <= 3; i ++){
        k_mer >>= 2;
        k_mer += (next[i] / 2) % 4;
        if(b_filter.contains(k_mer)){
            succs[i] = 1;
        }
    }
    return succs;
}

vector<bool> predecessor(kmer k_mer, BloomFilter b_filter){
    vector<bool> pred(4,0);
    string next = "ACGT";
    for(int i = 0; i <= 3; i ++){
        k_mer <<= 2;
        k_mer += (next[i] / 2) % 4;
        if(b_filter.contains(k_mer)){
            pred[i] = 1;
        }
    }
    return pred;
}

string build_forwards(string current_contig, kmer curr_kmer, BloomFilter index, int k){
    vector<bool> succs = successor(curr_kmer, index);
    vector<bool>::iterator it;
    if(accumulate(succs.begin(), succs.end(), 0) == 1){
        for(int i = 0; i <= 3; i++){
            if(succs[i] == 1){
                current_contig += "ACTG"[i];
                break;
            }
        }
        return build_forwards(current_contig, str2num(current_contig.substr(1,k)), index, k);
    }
    else{
        return current_contig;
    }
}

string build_backwards(string current_contig, kmer curr_kmer, BloomFilter index, int k){
    vector<bool> pred = predecessor(curr_kmer, index);
    vector<bool>::iterator it;
    if(accumulate(pred.begin(), pred.end(), 0) == 1){
        for(int i = 0; i <= 3; i++){
            if(pred[i] == 1){
                current_contig = "ACTG"[i] + current_contig;
                break;
            }
        }
        return build_backwards(current_contig, str2num(current_contig.substr(1,k)), index, k); 
    }
    else{
        return current_contig;
    }
}

string create_contig(string filename, int k, int h, int s){
    pair<BloomFilter, kmer> res = construct_index(filename, k, h, s);
    kmer curr_kmer = res.second;
    BloomFilter index = res.first;

    cout << "Starting first contig with the following " << to_string(k) << "-mer:" << kmer2str(curr_kmer, k) << endl;
    
    string contig = build_backwards(kmer2str(curr_kmer, k), curr_kmer, index, k)
    + build_forwards(kmer2str(curr_kmer, k), curr_kmer, index, k).substr(k);

    return contig;
}