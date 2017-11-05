#include <iostream>
#include <vector>
#include <fstream>
#include "MurmurHash3.h"
#include "BloomFilter.hpp"

#define NUM_HASH 200
#define LARGE_PRIME 9999999999971
#define FALSE_POSITIVE 0.001

using namespace std;

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter
vector<uint64_t> seeds;

vector<string> generate_shingles(string sequence, int k) {
    int size = sequence.size();
    vector <string> shingles;

    for (int i = 0; i <= size - k; i++) {
        shingles.push_back(sequence.substr(i, k));
    }

    return shingles;
}

double true_jaccard_similarity(string sequence1, string sequence2, int kmer_size) {
    int s1 = sequence1.size(), s2 = sequence2.size();
    vector <string> shingles1, shingles2, v;

    shingles1 = generate_shingles(sequence1, kmer_size);
    shingles2 = generate_shingles(sequence2, kmer_size);

    sort(shingles1.begin(), shingles1.end());
    sort(shingles2.begin(), shingles2.end());

    set_intersection(shingles1.begin(), shingles1.end(), shingles2.begin(), shingles2.end(), back_inserter(v));
    int num_intersection = v.size();

    v.clear();
    set_union(shingles1.begin(), shingles1.end(), shingles2.begin(), shingles2.end(), back_inserter(v));
    int num_union = v.size();

    return ((double)num_intersection) / num_union;
}

double min_hash_jaccard_estimate(vector <uint64_t> sketch1, vector <uint64_t> sketch2){
    int sketch_size = sketch1.size(), common = 0;

    for (int i = 0; i < sketch1.size(); i++) {
        if(sketch1[i] == sketch2[i]) {
            common++;
        }
    }

    return ((double)common / sketch_size);
}

double containment_jaccard_estimate(string sequence1, string sequence2, vector <string> sketch2, bloom_filter filter, int kmer_size) {
    int intersections = 0;
    int sketch_size = sketch2.size();
    for (int i = 0; i < sketch_size; i++) {
        if (filter.contains(sketch2[i])) {
            intersections++;
        }
    }
    intersections -= FALSE_POSITIVE * sketch_size;

    double containment_estimate = ((double)intersections / sketch_size);

    int size_sequence1_set = sequence1.size() - kmer_size + 1;
    int size_sequence2_set = sequence2.size() - kmer_size + 1;

    return ((double)(containment_estimate * size_sequence2_set)) / (size_sequence1_set + size_sequence2_set - size_sequence2_set * containment_estimate);
}

void generate_seeds() {
    srand (time(NULL));
    for (int i = 0; i < NUM_HASH; i++) {
        seeds.push_back(rand());
    }
}

uint64_t get_integer_fingerprint(string shingle, int hash_num) {
    const char *key = shingle.c_str();
    uint64_t hash_output[2];

    MurmurHash3_x64_128(key, (uint64_t)strlen(key), seeds[hash_num], hash_output);

    return (hash_output[0] +  NUM_HASH * hash_output[1]) % LARGE_PRIME;
}

vector<uint64_t> generate_sketch(vector <string> shingles) {
    int num_shingles = shingles.size();
    vector <uint64_t> sketch;

    for (int i = 0; i < NUM_HASH; i++) {
        uint64_t min_mer = LLONG_MAX;
        // cout << "Hash "<< i << ":\n---------------\n";
        for (int j = 0; j < num_shingles; j++) {
            uint64_t hash_value = get_integer_fingerprint(shingles[j], i);
            if (hash_value < min_mer) {
                min_mer = hash_value;
            }
            // cout << shingles[j] << " : " << hash_value << "\n";
        }
        // cout << "Min-mer : " << min_mer << "\n\n";
        sketch.push_back(min_mer);
    }

    return sketch;
}

vector<string> generate_kmer_sketch(vector <string> shingles) {
    int num_shingles = shingles.size();
    vector <string> sketch;

    for (int i = 0; i < NUM_HASH; i++) {
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

vector<string> read_dataset(string filename) {
    string sequence, line;
    vector <string> sequences;
    ifstream file (filename);

    if (file.is_open()) {
        getline(file, line);
        while (getline(file, line)) {
            if (line[0] != '>') {
                sequence += line;
            } else {
                sequences.push_back(sequence);
                sequence = "";
            }
        }
    } else {
        cout << "Unable to open file\n";
        exit(EXIT_FAILURE);
    }
    sequences.push_back(sequence);
    file.close();

    return sequences;
}

void print_sketch(vector <uint64_t> sketch) {
    int size = sketch.size();
    for (int i = 0; i < size; i++) {
        if (i == 0) {
            cout << "[ ";
        }
        cout << sketch[i];

        if (i == size - 1) {
            cout << " ]";
        } else {
            cout << ", ";
        }
    }
    cout << "\n";
}

double min_hash(string sequence1, string sequence2, int k) {
    double jaccard_index;
    vector <string> shingles1, shingles2;
    vector <uint64_t> sketch1, sketch2;

    shingles1 = generate_shingles(sequence1, k);
    sketch1   = generate_sketch(shingles1);

    shingles2 = generate_shingles(sequence2, k);
    sketch2   = generate_sketch(shingles2);

    // cout << "Sketch 1: ";
    // print_sketch(sketch1);
    // cout << "Sketch 2: ";
    // print_sketch(sketch2);

    return min_hash_jaccard_estimate(sketch1, sketch2);
}

double containment_hash(string sequence1, string sequence2, int kmer_size) {
    if (sequence2.size() > sequence1.size()) {
        string swap = sequence1;
        sequence1 = sequence2;
        sequence2 = swap;
    }

    bloom_parameters parameters;
    int s1 = sequence1.size();

    // How many elements roughly do we expect to insert?
    parameters.projected_element_count = s1 - kmer_size + 1; // Number of k mers

    // Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = FALSE_POSITIVE; // 1 in 1000

   // Simple randomizer (optional)
    parameters.random_seed = rand();

    if (!parameters) {
        cout << "Error - Invalid set of bloom filter parameters!" << endl;
        exit(EXIT_FAILURE);
    }

    parameters.compute_optimal_parameters();
    bloom_filter filter(parameters);

    string kmer;
    for (int i = 0; i <= s1 - kmer_size; i++) {
        kmer = sequence1.substr(i, kmer_size);
        if (!filter.contains(kmer)) {
            filter.insert(kmer);
        }
    }

    vector <string> shingles2;
    vector <string> sketch2;
    shingles2 = generate_shingles(sequence2, kmer_size);
    sketch2   = generate_kmer_sketch(shingles2);

    return containment_jaccard_estimate(sequence1, sequence2, sketch2, filter, kmer_size);
}

void sequences_similarity(vector <string> sequences1, vector <string> sequences2) {
    int n1 = sequences1.size(), n2 = sequences2.size(), kmer_size;
    double min_hash_jaccard, containment_hash_jaccard, jaccard_similarity;

    cout << "Enter value of k (Length of a shingle): ";
    cin >> kmer_size;
    generate_seeds();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            min_hash_jaccard = min_hash(sequences1[i], sequences2[j], kmer_size);
            containment_hash_jaccard = containment_hash(sequences1[i], sequences2[j], kmer_size);
            jaccard_similarity = true_jaccard_similarity(sequences1[i], sequences2[j], kmer_size);
            cout << "Min hash jaccard estimate: " << min_hash_jaccard << endl;
            cout << "Containment hash jaccard estimate: " << containment_hash_jaccard << endl;
            cout << "True jaccard similarity: " << jaccard_similarity << endl;
        }
    }
}

int main() {
    vector <string> sequences1, sequences2;
    sequences1 = read_dataset("dataset1.txt");
    sequences2 = read_dataset("dataset2.txt");

    sequences_similarity(sequences1, sequences2);

    return 0;
}
