#include <vector>

using namespace std;

#define LARGE_PRIME 9999999999971

void generate_seeds();

uint64_t get_integer_fingerprint(string shingle, int hash_num);

void generate_sketch(string shingle, vector<uint64_t> &min_sketch);

