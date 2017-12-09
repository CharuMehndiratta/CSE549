#include <vector>

using namespace std;

#define LARGE_PRIME 9999999999971

/* Number of hash functions, false positive and kmer size with default values */
double false_positive = 0.001;

int kmer_size = 16;

int num_hash = 100;

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter
vector<uint64_t> seeds;

void generate_seeds();

uint64_t get_integer_fingerprint(string shingle, int hash_num);

void generate_sketch(string shingle, vector<uint64_t> min_sketch);