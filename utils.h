#include <vector>
#include "bloomFilter.h"

using namespace std;

#define kmer uint64_t;

vector<string> find_kmers(const string& seq, int k);
kmer str2num(const string& str);
vector<bool> successor(kmer k_mer, BloomFilter b_filter);
vector<bool> predecessor(kmer k_mer, BloomFilter b_filter);
string kmer2str(kmer num, int l);