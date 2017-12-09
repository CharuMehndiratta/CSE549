#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include "MurmurHash3.h"
#include "BloomFilter.hpp"
#include <cstring>
#include <string.h>
#include <sstream>
#include <unistd.h>
#include "utils.h"

using namespace std;

/* Number of hash functions, false positive and kmer size with default values */
extern double false_positive;

extern int kmer_size;

extern int num_hash;

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter
extern vector<uint64_t> seeds;

string reference_file;

void read_min_sketch() {

    vector<uint64_t> min_sketch(num_hash, LLONG_MAX);
    ifstream file_min_read("min_sketch", ios::binary);

    file_min_read.read((char *)&min_sketch, sizeof(min_sketch));

    for(int i = 0; i < num_hash; i++){
        cout<<min_sketch[i]<<" ";
    }
    file_min_read.close();
}

void read_bloom_filter(string ref_genome) {

    bloom_filter ref_filter;
    ifstream file_bloom("bloom_filter", ios::binary);

    file_bloom.read((char *)&ref_filter, sizeof(ref_filter));
    int ref_size = ref_genome.size(), count = 0;

    for(int i = 0; i < ref_size - kmer_size + 1; i++){
        string kmer = ref_genome.substr(i, kmer_size);
            
        if (ref_filter.contains(kmer)) {
           count++;
        }
    }

    file_bloom.close();

    cout<<"\n count is "<<count;
}


void read_reference_genome(string reference_genome) {
    vector<uint64_t> reference_genome_min_sketch(num_hash, LLONG_MAX);
    string kmer;
    int ref_size = reference_genome.size();

    ofstream min_sketch_file("min_sketch", ios::binary);
    fstream bloom_filter_file("bloom_filter", ios::binary);

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
    int count = 0;

    for(int i = 0; i < ref_size - kmer_size; i++){
        kmer = reference_genome.substr(i, kmer_size);

        //Adding kmers to min hash
        generate_sketch(kmer, reference_genome_min_sketch);

        // Adding kmers to bloom filter
        if (!filter.contains(kmer)) {
            count++;
            filter.insert(kmer);
        }
    }

    min_sketch_file.write((char*)&reference_genome_min_sketch, sizeof(reference_genome_min_sketch));

    // Writing bloom filter to file
    bloom_filter_file.write((char*)&filter, sizeof(filter));

    read_bloom_filter(reference_genome);
    read_min_sketch();

}


void read_dataset(string filename) {

    string sequence, line;
    ifstream file (filename);

    if (file.is_open()) {
        getline(file, line);
        while (getline(file, line)) {
            if (line[0] != '>') {

                read_reference_genome(line);

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

    /* Get options from command line arguments */
    // while ((option = getopt(argc, argv, "r:f:k:h:")) != -1) {
    //     switch (option) {
    //         case 'r':
    //             reference_file = optarg;
    //             break;
    //         case 'f':
    //             false_positive = atof(optarg);
    //             break;
    //         case 'k':
    //             kmer_size = atoi(optarg);
    //             break;
    //         case 'h':
    //             num_hash = atoi(optarg);
    //             break;
    //         default:
    //             exit(EXIT_FAILURE);
    //     }
    // }

    reference_file = "reference_genome.txt";
    false_positive = 0.01;
    kmer_size = 16;
    num_hash = 10;

    generate_seeds();
    read_dataset(reference_file);


    return 0;
}
