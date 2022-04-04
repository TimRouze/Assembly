#include "utils.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <unistd.h>

class BloomFilter{
private:

    //Number of hash functions to apply
    size_t k;

    //Number of k-mers in the bloom filter
    size_t n = 0;

    //Size of the bloom filter
    size_t m;

    //The actual bloom filter
    vector<bool> filter;

public:

    //The xorshift hash function
    uint64_t xorshift(uint64_t x, uint64_t seed);

    //The constructor
    BloomFilter(size_t hash_func_number, size_t size);

    //To insert a k-mer in the filter
    void insert(uint64_t kmer);

    //To check if a k-mer is already in the filter
    bool contains(uint64_t kmer);

    //To get the number of k-mers in the filter
    int getN();
};

vector<string> filter(const vector<__uint128_t>& skmers, size_t k, size_t m, const vector<int>& sizes);
