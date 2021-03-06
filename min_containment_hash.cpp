#include <iostream>
#include <vector>
#include <fstream>
#include "MurmurHash3.h"
#include "BloomFilter.hpp"
#include <string.h>
#include <sstream>
#include <unistd.h>
#include <limits.h>

using namespace std;

#define LARGE_PRIME 9999999999971

/* Number of hash functions, false positive and kmer size with default values */
double false_positive = 0.001;

int kmer_size = 11;

int num_hash = 100;

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter
vector<uint64_t> seeds;

string reference_genome_file;

string reference_genome_min_sketch_file;

string reference_genome_bloom_filter_file;

string reference_genome_size_file;

string seeds_file;

string index_file;

string reference_genome_kmer_file;

void clear_file_content(string file) {
    ofstream file_to_clear(file, fstream::trunc);
    file_to_clear.close();
}

void clear_initial_files() {
    clear_file_content(index_file);
    clear_file_content(reference_genome_min_sketch_file);
    clear_file_content(reference_genome_bloom_filter_file);
    clear_file_content(reference_genome_size_file);
    clear_file_content(seeds_file);
    clear_file_content(reference_genome_kmer_file);
}

void set_initial_data() {
    reference_genome_min_sketch_file    = "reference_genome_min_sketch_file.txt";
    reference_genome_bloom_filter_file  = "reference_genome_bloom_filter_file";
    reference_genome_size_file          = "reference_genome_size_file.txt";
    reference_genome_kmer_file          = "reference_genome_kmer_file.txt";
    seeds_file                          = "seeds.txt";
    index_file                          = "index.txt";
}

void set_intermediate_file(){

    ofstream index_file_name;
    index_file_name.open(index_file);
    string line = "";

    //Write kmer_size on intermediate file
    line = "kmer_size " + to_string(kmer_size);
    index_file_name << line << "\n";
    line.clear();

    line = "false_positive " + to_string(false_positive);
    index_file_name << line << "\n";
    line.clear();

    line = "num_hash " + to_string(num_hash);
    index_file_name << line << "\n";
    line.clear();

}

void generate_seeds() {
    ofstream seed_file_ref;
    seed_file_ref.open(seeds_file);
    srand(time(NULL));
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

void read_bloom_filter(string ref_genome) {
    bloom_filter ref_filter;
    fstream file_bloom(reference_genome_bloom_filter_file, ios::binary);

    file_bloom.read((char *)&ref_filter, sizeof(ref_filter));
    int ref_size = ref_genome.size(), count = 0;

    for(int i = 0; i < ref_size - kmer_size + 1; i++){
        string kmer = ref_genome.substr(i, kmer_size);

        if (ref_filter.contains(kmer)) {
           count++;
        }
    }

    file_bloom.close();

}

void read_reference_genome(string reference_genome, ofstream &min_sketch_file, ofstream &bloom_filter_file, ofstream &reference_genome_kmer) {
    vector<uint64_t> reference_genome_min_sketch(num_hash, LLONG_MAX);
    string kmer;
    int ref_size = reference_genome.size();

    bloom_parameters parameters;

    // How many elements roughly do we expect to insert?
    parameters.projected_element_count = ref_size - kmer_size + 1; // Number of k mers

    // Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = false_positive; // 1 in 1000

   // Simple randomizer (optional)
    parameters.random_seed = rand();

    parameters.compute_optimal_parameters();

    if (!parameters) {
        exit(EXIT_FAILURE);
    }

    bloom_filter filter(parameters);

    for (int i = 0; i < ref_size - kmer_size; i++){
        kmer = reference_genome.substr(i, kmer_size);

        reference_genome_kmer << kmer <<" ";
        //Adding kmers to min hash
        generate_sketch(kmer, reference_genome_min_sketch);

        // Adding kmers to bloom filter
        if (!filter.contains(kmer)) {
            filter.insert(kmer);
        }
    }

    reference_genome_kmer << "\n";

    for(int i = 0; i < num_hash; i++) {
        min_sketch_file << reference_genome_min_sketch[i] << " ";
    }
    min_sketch_file << "\n";

    bloom_filter_file.write((char*)&filter, sizeof(filter));
}

void read_dataset(string filename) {
    ofstream ref_file(reference_genome_size_file);

    ofstream bloom_filter_file(reference_genome_bloom_filter_file, ios::binary | ios::app);
    ofstream min_sketch_file(reference_genome_min_sketch_file, ios::app);
    ofstream reference_genome_kmer(reference_genome_kmer_file, ios::app);

    ifstream file(filename);
    string sequence, line;

    string tmp = "";
    if (file.is_open()) {
        while (getline(file, tmp)) {
            if (tmp[0] == '>')
                continue;
            line = tmp;
            while (getline(file, tmp) && tmp[0] != '>')
                line += tmp;

            ref_file << to_string(line.size());
            ref_file << " ";

            read_reference_genome(line, min_sketch_file, bloom_filter_file, reference_genome_kmer);
            bloom_filter_file.flush();
            reference_genome_kmer.flush();
        }
    } else {
        cout << "Unable to open file\n";
        exit(EXIT_FAILURE);
    }
    ref_file.close();
    file.close();
    bloom_filter_file.close();
    min_sketch_file.close();
    reference_genome_kmer.close();
}

/*************************************************************/
/* Start execution                                           */
/*************************************************************/
int main(int argc, char *argv[]) {
    int option;

    /* Get options from command line arguments */
    while ((option = getopt(argc, argv, "r:f:k:h:")) != -1) {
        switch (option) {
            case 'r':
                reference_genome_file = optarg;
                break;
            case 'f':
                false_positive = atof(optarg);
                break;
            case 'k':
                kmer_size = atoi(optarg);
                break;
            case 'h':
                num_hash = atoi(optarg);
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }

    set_initial_data();
    clear_initial_files();
    generate_seeds();
    read_dataset(reference_genome_file);
    set_intermediate_file();

    return 0;
}
