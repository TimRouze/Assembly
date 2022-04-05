#include "utils.h"
#include "bloomFilter.h"
#include <numeric>
#include <string>
using namespace std;

BloomFilter construct_filter(string filename);
vector<bool> successor(kmer k_mer, BloomFilter b_filter);
vector<bool> predecessor(kmer k_mer, BloomFilter b_filter);
string build_forwards(string current_contig, kmer curr_kmer, BloomFilter index, int k);
string build_backwards(string current_contig, kmer curr_kmer, BloomFilter index, int k);
string create_contig(string filename, int k, int h, int s);