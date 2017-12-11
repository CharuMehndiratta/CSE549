#include <iostream>
#include "utils.h"
#include <vector>
#include <string.h>
#include <cstring>
#include "MurmurHash3.h"
#include <fstream>

using namespace std;

/* Number of hash functions, false positive and kmer size with default values */
double false_positive;

int kmer_size;

int num_hash;

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter
vector<uint64_t> seeds;

string reference_genome_file = "reference_genome.txt";

string reference_genome_min_sketch_file = "reference_genome_min_sketch_file.txt";

string reference_genome_bloom_filter_file = "reference_genome_bloom_filter_file";

string seeds_file = "seeds.txt";

string long_read_file = "long_read.txt";

string long_read_min_sketch_file = "long_read_min_sketch_file";

string min_hash_output = "min_hash_output.txt";

string containment_hash_output = "containment_hash_output.txt";

string long_read_containment_hash_file = "long_read_containment_hash_file.txt";

string reference_genome_size_file = "reference_genome_size_file.txt";

string seed_file = "seed_file.txt";


void generate_seeds() {

    ofstream seed_file_ref;
    seed_file_ref.open(seed_file);
    srand (time(NULL));
    for (int i = 0; i < num_hash; i++) {
        uint64_t val = rand();
        seeds.push_back(val);
        seed_file_ref << val << " ";
    }
    seed_file_ref << "\n";
    seed_file_ref.close();
}


uint64_t get_integer_fingerprint(string shingle, int hash_num) {
    const char *key = shingle.c_str();
    uint64_t hash_output[2];

    MurmurHash3_x64_128(key, (uint64_t)strlen(key), seeds[hash_num], hash_output);

    return (hash_output[0] +  num_hash * hash_output[1]) % LARGE_PRIME;

    return 1;

}

void generate_sketch(string shingle, vector<uint64_t> &min_sketch) {

    for (int i = 0; i < num_hash; i++) {
        uint64_t min_mer = LLONG_MAX;
        uint64_t hash_value = get_integer_fingerprint(shingle, i);

        if (hash_value < min_sketch[i]) {
            min_sketch[i] = hash_value;
        }
    }

}