#include "utils.h"

vector<string> find_kmers(const string& seq, int k){
    vector<string> kmers;
    for(int i = 0; i <= seq.length(); i++){
         if(seq.substr(i, k).length() == k){
             kmers.push_back(seq.substr(i, k));
         }
    }
    return kmers;
}
