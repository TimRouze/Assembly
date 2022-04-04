#include "bloomFilter.h"

//The xorshift hash function
uint64_t BloomFilter::xorshift(uint64_t x, uint64_t seed){
    uint64_t s = seed;
    x ^= x << (27 * s) ;
    x ^= x << (4 * s) ;
    x ^= x << (49 * s) ;
    return x;
}

BloomFilter::BloomFilter(size_t hash_func_number, size_t size){
    if(0 == hash_func_number)
    {
        throw invalid_argument("Bloomfilter could not be initialized: k must be larger than 0");
    }
    k = hash_func_number;
    m = size;
    filter = vector<bool>(size,false);
}

//To insert a k-mer in the filter
void BloomFilter::insert(uint64_t kmer){
    for (size_t i = 0; i < k; i++){
        int index = xorshift(kmer,i) % m;
        if (filter[index] == false){ 
            n++;
            filter[index] = true;
        }
    }
}

//To check if a k-mer is already in the filter
bool BloomFilter::contains(uint64_t kmer){
    for (size_t i = 0; i < k; i++){
        int index = xorshift(kmer, i) % m;
        if (!filter[index]) return false;
    }
    return true;
}

int BloomFilter::getN(){
    return n;
}

vector<string> filter(const vector<__uint128_t>& skmers, size_t k, size_t m, const vector<int>& sizes){
    string skmer;
    vector<string> filtered_skmers;
    BloomFilter* b_filter = new BloomFilter(k, m);
    for (int i = 0; i < skmers.size(); i++){
        if (!b_filter->contains(skmers[i])){
            b_filter->insert(skmers[i]);
            skmer = skmer2str(skmers[i], sizes[i]);
            filtered_skmers.push_back(skmer);
        }
    }
    cout << "Input file contained " << to_string(skmers.size()) << " kmers.\n";
    cout << "Bloom filter caught " << to_string(b_filter->getN()) << " unique kmers.\n";
    return filtered_skmers;
}