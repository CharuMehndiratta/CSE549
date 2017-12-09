#include <iostream>
#include "utils.h"
#include <vector>
#include <string.h>
#include <cstring>
#include "MurmurHash3.h"

using namespace std;

/* Number of hash functions, false positive and kmer size with default values */
double false_positive;

int kmer_size;

int num_hash;

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter
vector<uint64_t> seeds;

void generate_seeds() {
    srand (time(NULL));
    for (int i = 0; i < num_hash; i++) {
        seeds.push_back(rand());
    }
}

uint64_t get_integer_fingerprint(string shingle, int hash_num) {
    const char *key = shingle.c_str();
    uint64_t hash_output[2];

    MurmurHash3_x64_128(key, (uint64_t)strlen(key), seeds[hash_num], hash_output);

    return (hash_output[0] +  num_hash * hash_output[1]) % LARGE_PRIME;
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