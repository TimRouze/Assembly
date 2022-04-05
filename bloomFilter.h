#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <stdexcept>
#include <vector>

using namespace std;

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

    //Size of the block window
    const size_t w = 2^8;

public:

    //The xorshift hash function
    uint64_t xorshift(uint64_t x, uint64_t seed);

    //The constructor
    BloomFilter(size_t hash_func_number, size_t size);

    void blocked_insert(uint64_t kmer);

    //To insert a k-mer in the filter
    void insert(uint64_t kmer);

    bool blocked_contains(uint64_t kmer);

    //To check if a k-mer is already in the filter
    bool contains(uint64_t kmer);

    //To get the number of k-mers in the filter
    int getN();
};
