int main(int argc, char** argv){
    char ch;
    string filename, output_filename;
	int r = 10;
    int k = 31;
    int m = 8;
	int h = 2;
	int s = 1000000;
    while ((ch = getopt(argc, argv, "q:k:f:m:h:s:r:")) != -1) {
		switch (ch) {
			case 'q': filename = optarg; break;
			case 'f': output_filename = optarg; break;
			case 'r': r = stoi(optarg); break;
			case 'k': k = stoi(optarg); break;
            case 'm': m = stoi(optarg); break;
			case 'h': h = stoi(optarg); break;
			case 's': s = stoi(optarg); break;
		}
	}
    if (output_filename == "" && filename != "") {
		output_filename = "results.fa";
	}
	if ((filename == "" && output_filename == "") | k == 0) {
		cout << "Core arguments:" << endl
		     << "	-q query file" << endl
             << "	-f output file" << endl
			 << "	-r subsampling rate" << endl
		     << "	-k k value used for kmers (31) " << endl
		     << "	-m minimizer size (8)" << endl
			 << "	-h number of hash function used for the bloom filter (2)" << endl
			 << "	-s bloom filter size (1Go)" << endl;
		return 0;
	}

    {
		cout << "Values used -k " + to_string(k) + " -r" + to_string(r) + " -m  " + to_string(m) + " -f " + 
				output_filename + " -h " + to_string(h) + " -s " + to_string(s)
		     << endl;
		if (filename != "" && output_filename != "") {
			subsampler(filename, output_filename, r, k, m, h, s);
		}
    }
    

    return 0;
