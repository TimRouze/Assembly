#include "utils.h"
#include "bloomFilter.h"
#include <numeric>
#include <string>
using namespace std;

pair <BloomFilter, string> construct_index(string filename, int k, int h, int s);
vector<bool> successor(string k_mer, BloomFilter b_filter);
vector<bool> predecessor(string k_mer, BloomFilter b_filter);
string build_forwards(string current_contig, string curr_kmer, BloomFilter index, int cpt);
string build_backwards(string current_contig, string curr_kmer, BloomFilter index);
string create_contig(string filename, int k, int h, int s);