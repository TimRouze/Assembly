#include "assembler.h"

BloomFilter read_file(string filename, int k, int s){
    ifstream fasta(filename);
    BloomFilter * b_filter = new BloolFilter(k, s);
    string sequence = "";
	if(fasta.good()){
		string desc;
		while(not fasta.eof()){
			getline(fasta, desc);
			getline(fasta, line);
            for(int i = 0; i < line.length(); i++){
                string k_mer = 
            }
        }
    }
    return sequence;
}

void assemble(string filename, int k, int h, int s){
    string sequence = read_file(filename, k)
}