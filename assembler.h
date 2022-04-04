#include "utils.h"
#include "bloomFilter.h"
using namespace std;

BloomFilter construct_filter(string filename);
void create_contig(string filename, int k, int h, int s);