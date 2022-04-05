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

void BloomFilter::insert(uint64_t kmer){
    bool first_insertion;
    int first_index = xorshift(kmer,1) % (m - w);
    first_insertion = !filter[first_index];
    filter[first_index] = true;
    for (size_t i = 2; i <= k; i++){
        int index = first_index + xorshift(kmer,i) % w;
        first_insertion |= !filter[index];
        filter[index] = true;
    }
    if (first_insertion){
        n++;
    }
}

bool BloomFilter::contains(uint64_t kmer){
    int first_index = xorshift(kmer,1) % (m - w);
    if (!filter[first_index]) return false;
    for (size_t i = 2; i <= k; i++){
        int index = first_index + xorshift(kmer,i) % w;
        if (!filter[index]) return false;
    }
    return true;
}

int BloomFilter::getN(){
    return n;
}