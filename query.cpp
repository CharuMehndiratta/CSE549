#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <sstream>
#include "MurmurHash3.h"
#include "utils.h"
#include "BloomFilter.hpp"

using namespace std;

/* Number of hash functions, false positive and kmer size with default values */
extern double false_positive;

extern int kmer_size;

extern int num_hash;

extern string min_hash_output;

extern string reference_genome_min_sketch_file;

extern string reference_genome_bloom_filter_file;

extern string reference_genome_size_file;

extern string long_read_file;

extern string containment_hash_output;

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

    ofstream min_hash_output_file(min_hash_output);

    // Reading reference genome min hash sketches 
    ifstream genome_sketch(reference_genome_min_sketch_file, ios::binary);

    while(genome_sketch.read((char *)&genome_min_sketch, sizeof(genome_min_sketch)) != NULL) {
        min_hash_jaccard_index =  min_hash_jaccard_estimate(genome_min_sketch, long_read_sketch);
        min_hash_output_file << min_hash_jaccard_index << " ";
    }
}

double containment_jaccard_estimate(int sequence1_size, string sequence2, vector <string> sketch2, bloom_filter filter) {
    int intersections = 0;
    int sketch_size = sketch2.size();
    for (int i = 0; i < sketch_size; i++) {
        if (filter.contains(sketch2[i])) {
            intersections++;
        }
    }
    intersections -= false_positive * sketch_size;

    double containment_estimate = ((double)intersections / sketch_size);

    int size_sequence1_set = sequence1_size - kmer_size + 1;
    int size_sequence2_set = sequence2.size() - kmer_size + 1;

    return ((double)(containment_estimate * size_sequence2_set)) / (size_sequence1_set + size_sequence2_set - size_sequence2_set * containment_estimate);
}

vector<string> generate_kmer_sketch(vector <string> shingles) {
    int num_shingles = shingles.size();
    vector <string> sketch;

    for (int i = 0; i < num_hash; i++) {
        uint64_t min_mer = LLONG_MAX;
        string kmer;
        for (int j = 0; j < num_shingles; j++) {
            uint64_t hash_value = get_integer_fingerprint(shingles[j], i);
            if (hash_value < min_mer) {
                min_mer = hash_value;
                kmer = shingles[j];
            }
        }
        sketch.push_back(kmer);
    }

    return sketch;
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
    bloom_filter ref_bloom_filter;
    int count = 0;
    ifstream ref_genome_size_file(reference_genome_size_file);
    string line, kmer;
    double jaccard_estimate;
    ofstream containment_hash_output_file(containment_hash_output);

    getline(ref_genome_size_file, line);

    ifstream bloom_filter_file(reference_genome_bloom_filter_file, ios::binary);

    bloom_filter_file.read((char *)&ref_bloom_filter, sizeof(ref_bloom_filter));

    // generate long read shingles
    vector <string> shingles_long_read;
    vector <string> sketch_long_read;

    shingles_long_read = generate_shingles(long_read, kmer_size);
    sketch_long_read   = generate_kmer_sketch(shingles_long_read);

    stringstream lines(line);

    int size;

    while (lines >> size) {
        jaccard_estimate = containment_jaccard_estimate(size, long_read, sketch_long_read, ref_bloom_filter);
        containment_hash_output_file << jaccard_estimate << " ";
    }

}

void generate_jacard_index(string long_read) {


    //true_jacard(long_read);
    // min_hash(long_read);
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
    // char *long_read_file = NULL, *bloom_filter_file = NULL, *min_sketch_file = NULL;

    /* Get options from command line arguments */
    // while ((option = getopt(argc, argv, "r:b:m")) != -1) {
    //     switch (option) {
    //         case 'r':
    //             long_read_file = optarg;
    //             break;
    //         case 'b':
    //             bloom_filter_file = optarg;
    //             break;
    //         case 'm':
    //             min_sketch_file = optarg;
    //             break;
    //         default:
    //             exit(EXIT_FAILURE);
    //     }
    // }

    read_dataset(long_read_file);

    return 0;
}
