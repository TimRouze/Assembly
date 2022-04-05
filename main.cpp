#include "assembler.h"

int main(int argc, char** argv){
    char ch;
    string filename;
    int k = 31;
	int h = 4;
	int s = 5000000;
    while ((ch = getopt(argc, argv, "q:k:h:s")) != -1) {
		switch (ch) {
			case 'q': filename = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 'h': h = stoi(optarg); break;
			case 's': s = stoi(optarg); break;
		}
	}
	if (filename == "" || k == 0 || h == 0){
		cout << "Core arguments:" << endl
		     << "	-q sequencing file" << endl
		     << "	-k kmers word size (31) " << endl
			 << "	-h number of hash function used for the bloom filter (4)" << endl
			 << "	-s bloom filter size (1Gb)" << endl;
		return 0;
	}
	else{
		cout << "Values used -q " + filename + " -k " + to_string(k) + " -h " + to_string(h) + " -s " + to_string(s)
		     << endl;
		if (filename != "") {
			string contig = create_contig(filename, k, h, s);
			cout << contig << endl;
			cout << contig.length() << endl;
		}
    }
    

    return 0;
}
