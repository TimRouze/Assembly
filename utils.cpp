#include "utils.h"

vector<string> find_kmers(const string& seq, int k){
    vector<string> kmers;
    for(int i = 0; i <= seq.length(); i++){
         if(seq.substr(i, k).length() == k){
             kmers.push_back(seq.substr(i, k));
         }
    }
    return kmers;
}

kmer str2num(const string& str) {
	kmer res(0);
	for (uint64_t i(0); i < str.size(); i++) {
		res <<= 2;
		res += (str[i] / 2) % 4;
	}
	return res;
}

vector<bool> successor(kmer k_mer, BloomFilter b_filter){
    vector<bool> succs = [0, 0, 0, 0];
    string next = "ACGT";
    for(int i = 0; i <= 3; i ++){
        k_mer >> 2;
        k_mer += (str[i] / 2) % 4;
        if(b_filter -> contains(k_mer)){
            succs[i] = 1;
        }
    }
    return succs;
}

vector<bool> predecessor(kmer k_mer, BloomFilter b_filter){
    vector<bool> preds = [0, 0, 0, 0];
    string next = "ACGT";
    for(int i = 0; i <= 3; i ++){
        k_mer << 2;
        k_mer += (str[i] / 2) % 4;
        if(b_filter -> contains(k_mer)){
            preds[i] = 1;
        }
    }
    return preds;
}

string kmer2str(kmer num, int l){
	string res(l, '\0');
	Pow2<kmer> anc(2 * (l - 1));
	for (uint64_t i(0); i < l; ++i) {
		uint64_t nuc = num / anc;
		num          = num % anc;
		assert(nuc < 4);
		res[i] = "ACTG"[nuc];
		anc >>= 2;
	}
	return res;
}

