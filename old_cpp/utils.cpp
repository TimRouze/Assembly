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

