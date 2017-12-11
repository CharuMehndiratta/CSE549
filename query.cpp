#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <sstream>
#include "MurmurHash3.h"
// #include "utils.h"
#include "BloomFilter.hpp"

using namespace std;

/* Number of hash functions, false positive and kmer size with default values */
 double false_positive;

 int kmer_size;

 int num_hash;

vector<uint64_t> seeds;

string min_hash_output;

string reference_genome_min_sketch_file;

string reference_genome_bloom_filter_file;

string reference_genome_size_file;

string long_read_file;

string containment_hash_output;

string seeds_file;

string index_file;

#define LARGE_PRIME 9999999999971

// https://stackoverflow.com/questions/9241230/what-is-murmurhash3-seed-parameter




void clear_file_content(string file){
    ofstream file_to_clear(file, fstream::trunc);
    file_to_clear.close();
}




void clear_initial_files(){
    clear_file_content(min_hash_output);
    clear_file_content(containment_hash_output);
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


void read_index_file(){
    ifstream ref_index_file(index_file);
    string line = "";
    while(getline(ref_index_file, line)){
        stringstream ss(line); 
        string val;
        ss >> val;
        if(val.compare("kmer_size") == 0){
            ss >> kmer_size;
        }
        else if(val.compare("false_positive") == 0){
            ss >> false_positive;
        }
        else if(val.compare("num_hash") == 0){
            ss >> num_hash;
        }
        else {
            cout << "index_file corrupted\n";
        }

    }


}



void read_seed_file() {

    ifstream seed_file_ref(seeds_file);
    for(int i = 0; i < num_hash; i++ ) {
        uint64_t val;
        seed_file_ref >> val;
        seeds.push_back(val);
    }
    seed_file_ref.close();
}

/*************************************************************/
/*  Input:   long_read and reference hash signature          */
/*  output:   jaccard index between the input signatures.    */
/*************************************************************/

double min_hash_jaccard_estimate(vector <uint64_t> sketch1, vector <uint64_t> sketch2){

    int sketch_size = sketch1.size(), common = 0;

    for (int i = 0; i < sketch1.size(); i++) {
        if(sketch1[i] == sketch2[i]) {
            common++;
        }
    }

    return ((double)common / sketch_size);
}





void set_initial_data(){

    reference_genome_min_sketch_file = "reference_genome_min_sketch_file.txt";
    reference_genome_bloom_filter_file  = "reference_genome_bloom_filter_file";
    reference_genome_size_file  = "reference_genome_size_file.txt";
    seeds_file  = "seeds.txt";
    min_hash_output = "min_hash_output.txt";
    containment_hash_output = "containment_hash_output.txt";
}



/*************************************************************/
/*  Input:   min_hash signature reference, long_read         */
/*  output:  min_hash signature on all k-mers of long_read   */
/*************************************************************/
void gen_read_sketch(vector<uint64_t> & sketch, string read){

    int size = read.size();
    for(int i = 0; i < size - kmer_size + 1; i++){
        string kmer = read.substr(i, kmer_size);
        generate_sketch(kmer, sketch);
    }

    cout << "\n\n";
}


/*************************************************************/
/*  Input: long read                                         */
/*  output: jaccard index of long_read input with all the.   */
/*    genome references                                      */
/*************************************************************/
void min_hash(string long_read){
        
    vector<uint64_t> long_read_sketch(num_hash, LLONG_MAX);
    vector<uint64_t> genome_min_sketch;
    double min_hash_jaccard_index;
    gen_read_sketch(long_read_sketch, long_read);
    
    string line = "";
    
    fstream min_hash_output_file(min_hash_output, ios::app);

    ifstream genome_sketch(reference_genome_min_sketch_file, ios::binary);

    while(getline(genome_sketch, line)){

        stringstream liness(line);
        uint64_t val;
        while(liness >> val){
            genome_min_sketch.push_back(val);
        }

        min_hash_jaccard_index = min_hash_jaccard_estimate(long_read_sketch, genome_min_sketch);

        min_hash_output_file <<  min_hash_jaccard_index << " ";
        genome_min_sketch.clear();
    }
    min_hash_output_file << "\n";
    min_hash_output_file.close();
    genome_sketch.close();
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

/*************************************************************/
/*  Input: long read, kmer size                              */
/*  output: kmer sized shingles stored in a vector<string>   */
/*************************************************************/
vector<string> generate_shingles(string sequence) {
    int size = sequence.size();
    vector <string> shingles;
    string kmer;

    for (int i = 0; i <= size - kmer_size; i++) {
        kmer = sequence.substr(i, kmer_size);
        shingles.push_back(kmer);
    }

    return shingles;
}




void containment_hash(string long_read) {

    bloom_filter ref_bloom_filter;
    int count = 0, size;
    string line, kmer;
    double jaccard_estimate;

    ifstream ref_genome_size_file(reference_genome_size_file);
    ofstream containment_hash_output_file(containment_hash_output, ios::app);
    fstream bloom_filter_file(reference_genome_bloom_filter_file, ios::binary);

    getline(ref_genome_size_file, line);

    bloom_filter_file.read((char *)&ref_bloom_filter, sizeof(ref_bloom_filter));

    // generate long read shingles
    vector <string> shingles_long_read;
    vector <string> sketch_long_read;

    shingles_long_read = generate_shingles(long_read);
    sketch_long_read   = generate_kmer_sketch(shingles_long_read);

    stringstream lines(line);

    while (lines >> size) {
        jaccard_estimate = containment_jaccard_estimate(size, long_read, sketch_long_read, ref_bloom_filter);
        containment_hash_output_file << to_string(jaccard_estimate) << " ";
    }
    containment_hash_output_file << "\n";
    bloom_filter_file.close();
    ref_genome_size_file.close();
    containment_hash_output_file.close();

}



/*************************************************************/
/*   This function takes a long read as input and calls      */
/*    min_hash and containment hash function to calculate    */
/*    jaccard similarity                                     */
/*************************************************************/
void generate_jacard_index(string long_read) {
    //true_jacard(long_read);
    min_hash(long_read);
    containment_hash(long_read);
}


void read_dataset(string filename) {

    string sequence, line, tmp = "";
    ifstream file (filename);

    if (file.is_open()) {
        while (getline(file, tmp)) {
            if (tmp[0] == '>')
                continue;

            line = tmp;
            while(getline(file, tmp) && tmp[0] != '>')
                line += tmp;

            generate_jacard_index(line);
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
   
     while ((option = getopt(argc, argv, "r:i:")) != -1) {
        switch (option) {
            case 'r':
                long_read_file = optarg;
                break;
            case 'i':
                index_file = optarg;
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }

    set_initial_data();
    clear_initial_files();
    read_index_file();
    read_seed_file();

    read_dataset(long_read_file);

    return 0;
}
