#include "assembler.h"

pair <BloomFilter, string> construct_index(string filename, int k, int h, int s){
    ifstream fasta(filename);
    BloomFilter *index = new BloomFilter(h, s);
    string sequence = "";
    string k_mer;
    int n = 0;
    if (not fasta.good()) {
		cout << "Problem with file opening" << endl;
		exit(1);
	}
	else{
		string desc;
		while(not fasta.eof()){
			getline(fasta, desc);
			getline(fasta, sequence);
            for(int i = 0; i < sequence.length() - k + 1; i++){
                k_mer = sequence.substr(i, k);
                n++;
                if (!index->contains(str2num(k_mer))){
                    index->insert(str2num(k_mer));
                }
            }
        }
    }
    string file = filename.substr(filename.find_last_of("/")+1, filename.length());
    cout << file + " contains " + to_string(n) + " " + to_string(k) + "-mers."  << endl;
    cout << "The Bloom Filter caught " + to_string(index->getN()) + " unique " + to_string(k) + "-mers." << endl;

    pair<BloomFilter, string> res (*index, k_mer);
    return res;
}

vector<bool> successor(string k_mer, BloomFilter b_filter){
    vector<bool> succs(4,0);
    string next = "ACTG";
    k_mer.erase(k_mer.begin());
    for(int i = 0; i <= 3; i ++){
        k_mer += next[i];
        /*cout << k_mer << endl;
        cout << b_filter.contains(str2num(k_mer)) << endl;
        cin.get();*/
        if(b_filter.contains(str2num(k_mer))){
            succs[i] = 1;
        }
        k_mer.pop_back();
    }
    return succs;
}

vector<bool> predecessor(string k_mer, BloomFilter b_filter){
    vector<bool> pred(4,0);
    string next = "ACTG";
    k_mer.pop_back();
    for(int i = 0; i <= 3; i ++){
        k_mer = next[i] + k_mer;
        if(b_filter.contains(str2num(k_mer))){
            pred[i] = 1;
        }
        k_mer.erase(k_mer.begin());
    }
    return pred;
}

string build_forwards(string current_contig, string curr_kmer, BloomFilter index, int cpt){
    vector<bool> succs = successor(curr_kmer, index);
    if(accumulate(succs.begin(), succs.end(), 0) == 1){
        for(int i = 0; i <= 3; i++){
            if(succs[i] == 1){
                current_contig += "ACTG"[i];
                cpt += 1;
                /*cout << cpt << endl;
                cout << "next k_mer: " << current_contig.substr(cpt,curr_kmer.length()+1) << endl;
                cin.get();*/
                break;
            }
        }
        return build_forwards(current_contig, current_contig.substr(cpt,curr_kmer.length()+1), index, cpt);
    }
    else{
        cout << "kmer a traiter: " << curr_kmer << endl;
        cout << succs[0] << succs[1] << succs[2] << succs[3] << endl;
        cin.get();
        return current_contig;
    }
}

string build_backwards(string current_contig, string curr_kmer, BloomFilter index){
    vector<bool> pred = predecessor(curr_kmer, index);
    if(accumulate(pred.begin(), pred.end(), 0) == 1){
        for(int i = 0; i <= 3; i++){
            if(pred[i] == 1){
                current_contig = "ACTG"[i] + current_contig;
                break;
            }
        }
        return build_backwards(current_contig, current_contig.substr(0, curr_kmer.length()), index); 
    }
    else{
        return current_contig;
    }
}

string create_contig(string filename, int k, int h, int s){
    pair<BloomFilter, string> res = construct_index(filename, k, h, s);
    string curr_kmer = res.second;
    BloomFilter index = res.first;

    cout << "Starting first contig with the following " << to_string(k) << "-mer: " << curr_kmer << endl;

    string contig = build_backwards(curr_kmer, curr_kmer, index)
    + build_forwards(curr_kmer, curr_kmer, index, 0).substr(k);

    return contig;
}