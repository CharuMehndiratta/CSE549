#include <iostream>
#include <vector>
#include <fstream>
#include "MurmurHash3.h"

#define NUM_HASH 200
#define MAX_HASH_SIZE 1000000000

using namespace std;

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter
vector<uint64_t> seeds;

double jaccard_similarity(vector <uint64_t> sketch1, vector <uint64_t> sketch2){
    int sketch_size = sketch1.size(), common = 0;

    for (int i = 0; i < sketch1.size(); i++) {
        if(sketch1[i] == sketch2[i]) {
            common++;
        }
    }

    return ((double)common / sketch_size);
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

    return (hash_output[0] +  NUM_HASH * hash_output[1]) % MAX_HASH_SIZE;
}

vector<string> generate_shingles(string sequence, int k) {
    int size = sequence.size();
    vector <string> shingles;

    for (int i = 0; i <= size - k; i++) {
        shingles.push_back(sequence.substr(i, k));
    }

    return shingles;
}

vector<uint64_t> generate_sketch(vector <string> shingles) {
    int num_shingles = shingles.size();
    vector <uint64_t> sketch;

    for (int i = 0; i < NUM_HASH; i++) {
        uint64_t min_mer = LLONG_MAX;
        cout << "Hash "<< i << ":\n---------------\n";
        for (int j = 0; j < num_shingles; j++) {
            uint64_t hash_value = get_integer_fingerprint(shingles[j], i);
            if (hash_value < min_mer) {
                min_mer = hash_value;
            }
            cout << shingles[j] << " : " << hash_value << "\n";
        }
        cout << "Min-mer : " << min_mer << "\n\n";
        sketch.push_back(min_mer);
    }

    return sketch;
}

vector<string> read_dataset() {
    string sequence, line;
    vector <string> sequences;
    ifstream file ("dataset.txt");

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

void min_hash(string sequence1, string sequence2) {
    int k;
    double jaccard_index;
    vector <string> shingles1, shingles2;
    vector <uint64_t> sketch1, sketch2;

    cout << "Enter value of k: ";
    cin >> k;

    generate_seeds();
    shingles1 = generate_shingles(sequence1, k);
    sketch1   = generate_sketch(shingles1);

    shingles2 = generate_shingles(sequence2, k);
    sketch2   = generate_sketch(shingles2);

    cout << "Sketch 1: ";
    print_sketch(sketch1);
    cout << "Sketch 2: ";
    print_sketch(sketch2);

    jaccard_index = jaccard_similarity(sketch1, sketch2);

    cout << "Jaccard index: " << jaccard_index;
}

void sequence_similarity(vector <string> sequences) {
    int num_sequences = sequences.size();

    for (int i = 0; i < num_sequences - 1; i++) {
        for (int j = 1; j < num_sequences; j++) {
            min_hash(sequences[i], sequences[j]);
        }
    }
}

int main() {
    vector <string> sequences;
    sequences = read_dataset();

    sequence_similarity(sequences);

    return 0;
}
