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

void BloomFilter::blocked_insert(uint64_t kmer){
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

bool BloomFilter::blocked_contains(uint64_t kmer){
    int first_index = xorshift(kmer,1) % (m - w);
    if (!filter[first_index]) return false;
    for (size_t i = 2; i <= k; i++){
        int index = first_index + xorshift(kmer,i) % w;
        if (!filter[index]) return false;
    }
    return true;
}

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

 /*int main(){
    ifstream fasta("Data/Sequencage1_sans_erreurs.fa");
    BloomFilter *b = new BloomFilter(10,1000000);
    string sequence = "";
    int n_total_kmers = 0;
    if(fasta.good()){
        string desc;
        while(not fasta.eof()){
            getline(fasta, desc);
            getline(fasta, sequence);
            for(int i=0; i<sequence.length()-31+1; i++){
                kmer k_mer = str2num(sequence.substr(i, 31));
                n_total_kmers+=1;
                if(!b->blocked_contains(k_mer)){
                    b->blocked_insert(k_mer);
                }
            }
        }
    }
    cout << "Input file contained " << to_string(n_total_kmers) << " kmers.\n";
    cout << "Bloom filter caught " << to_string(b->getN()) << " unique kmers.\n";
}*/