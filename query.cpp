#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include "MurmurHash3.h"

using namespace std;

/* Number of hash functions, false positive and kmer size with default values */
extern double false_positive;

extern int kmer_size;

extern int num_hash;

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter
extern vector<uint64_t> seeds;


double min_hash_jaccard_estimate(vector <uint64_t> sketch1, vector <uint64_t> sketch2){
    int sketch_size = sketch1.size(), common = 0;

    for (int i = 0; i < sketch1.size(); i++) {
        if(sketch1[i] == sketch2[i]) {
            common++;
        }
    }

    return ((double)common / sketch_size);
}

void gen_read_sketch(vector<uint64_t> & sketch, string read){

    int size = read.size();
    for(int i = 0; i < size - num_hash + 1; i++){
        string kmer = read.substr(i, kmer_size);
        generate_sketch(kmer, sketch);
    }
}

// Input long read and generates jaccard index using min hash
void min_hash(string long_read){
    vector<uint64_t> long_read_sketch(num_hash, LLONG_MAX);
    vector<uint64_t> genome_min_sketch;
    double min_hash_jaccard_index;
    gen_read_sketch(long_read_sketch, long_read);

    ofstream min_hash(min_hash_output);

    // Reading reference genome min hash sketches 
    ifstream genome_sketch(reference_genome_min_sketch_file, ios::binary);

    while(genome_sketch.read((char *)&genome_min_sketch, sizeof(genome_min_sketch)) > 0) {
        min_hash_jaccard_index =  min_hash_jaccard_estimate(genome_min_sketch, long_read_sketch);
        min_hash << min_hash_jaccard_index << " ";
    }
}

vector<string> generate_shingles(string sequence, int kmer_size) {
    int size = sequence.size();
    vector <string> shingles;

    for (int i = 0; i <= size - kmer_size; i++) {
        shingles.push_back(sequence.substr(i, kmer_size));
    }

    return shingles;
}

void containment_hash(string long_read) { 
}

void generate_jacard_index(string long_read) {


    //true_jacard(long_read);
    min_hash(long_read);
    containment_hash(long_read);

}

void read_dataset(string filename) {

    string sequence, line;
    ifstream file (filename);

    if (file.is_open()) {
        getline(file, line);
        while (getline(file, line)) {
            if (line[0] != '>') {

                generate_jacard_index(line);

            //     sequence += line;
            // } else {
            //     cout<<"\n hello";
                
            //     sequence = "";
            }
        }
    } else {
        cout << "Unable to open file\n";
        exit(EXIT_FAILURE);
    }
    file.close();
}


/*************************************************************/
/* Start execution                                           */
/*************************************************************/
int main(int argc, char *argv[]) {
    int option;
    char *long_read_file = NULL, *bloom_filter_file = NULL, *min_sketch_file = NULL;

    /* Get options from command line arguments */
    while ((option = getopt(argc, argv, "r:b:m")) != -1) {
        switch (option) {
            case 'r':
                long_read_file = optarg;
                break;
            case 'b':
                bloom_filter_file = optarg;
                break;
            case 'm':
                min_sketch_file = optarg;
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }

    if (long_read_file == NULL) {

    }

    read_dataset(long_read_file);

    return 0;
}
