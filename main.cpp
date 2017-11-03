#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <string>
#include <cmath>
#include "MurmurHash3.h"
#include "BloomFilter.hpp"

#define MAX_HASH_ALGO 200
#define MAX_HASH_SIZE 10000000000
#define FALSE_POSTIVITY 0.001

using namespace std;

float jaccard_similarity(vector<uint64_t> & read1, vector<uint64_t>& read2){
    float common = 0, total = read1.size();
    float jaccard;

    for (int i = 0; i < read1.size(); i++) {
        if(read1[i] == read2[i]) {
            common++;
        }
    }

    cout << common << " / " << total <<"\n";
    jaccard = common/total;
    return jaccard;
}


vector<uint64_t> random_hash_vector;
void generate_rest_random_hash() {
    srand (time(NULL));
    for (int i = 0; i < MAX_HASH_ALGO; i++) {
        int x = rand();
        cout<<"\n rand -- "<<i;
        random_hash_vector.push_back(i);
    }

}


void generate_hash_value(int k, string str, vector<uint64_t> &hash_value, vector<string> &hash_string) {
    
    uint64_t hash_otpt[2], temp_hash;
    string temp;
    cout<<"\n input  "<<str;
    cout<<"\n";
    for (int i = 0; i <= str.size() - k; i++) {
        temp = str.substr(i, k);
        const char *key = temp.c_str();
        cout<<"\n key - "<<temp;
        for( int j = 0; j < MAX_HASH_ALGO; j++) {
            MurmurHash3_x64_128(key, (uint64_t)strlen(key), random_hash_vector[j], hash_otpt);
            temp_hash = (hash_otpt[0] +  MAX_HASH_ALGO * hash_otpt[1]) % MAX_HASH_SIZE;
            cout<<" temp_hash "<<temp_hash;

            if(hash_value[j] > temp_hash) {
                hash_value[j] = temp_hash;
                hash_string[j] = temp;
            }
        }
    }
}


int main() {
    srand (time(NULL));
    string input = "", tmp;
    vector < vector<uint64_t> > hash_values;
    vector < vector<string> > hash_strings;
    vector < string > strings;
    int n;
    ifstream file ("dataset.txt");
    int kmerSize = 3;
    
    if (file.is_open()) {
        getline(file, tmp);
        while(getline(file, tmp)) {
            if(tmp[0] != '>') {
                input += tmp;
            }
            else {
                strings.push_back(input);
                input = "";
            }
        }
    } else {
        cout<<"\n Unable to open file";
    }
    file.close();
    
    strings.push_back(input);
    
    n = int(strings.size());
    generate_rest_random_hash();
    for(int i = 0; i < n; i++) {
        vector<uint64_t> tmp(MAX_HASH_ALGO, LLONG_MAX);
        vector<string> temp_str(MAX_HASH_ALGO);
        hash_values.push_back(tmp);
        hash_strings.push_back(temp_str);
        generate_hash_value(kmerSize, strings[i], hash_values[i], hash_strings[i]);
        
        cout<<"\n-------------\n";
        for(int j = 0; j < MAX_HASH_ALGO; j++) {
            cout<<hash_values[i][j]<<" => ";
            cout<<hash_strings[i][j]<<" ";
        }
    }
    
    for(int k = 0; k < hash_values.size()-1; k++){
        for(int l = k + 1; l < hash_values.size(); l++){
            cout << "\nJaccard similarity between index " << k << " and index " << l << " using min hash approach is: ";
            cout << jaccard_similarity(hash_values[k], hash_values[l]);
            cout << "\n";
        }
    }
    
    int capacity = 2 * strings[0].size();

    bloom_parameters parameters;

    // How many elements roughly do we expect to insert?
    parameters.projected_element_count = capacity;

    // Maximum tolerable false positive probability? (0,1)
    parameters.false_positive_probability = FALSE_POSTIVITY ;

    // Simple randomizer (optional)
    parameters.random_seed = rand();
    cout<<"\n parameters.random_seed  "<<parameters.random_seed ;

    if (!parameters) {
      cout << "Error - Invalid set of bloom filter parameters!";
      return 1;
    }
    parameters.compute_optimal_parameters();
    bloom_filter filter(parameters);
    string temp;
    int size_B_est = 0; // String 0
    int size_A = 0; // String 1

    //Inserting String 0 in Bloom Filter
    for (int i = 0; i <= strings[0].size() - kmerSize; i++) {
        temp = strings[0].substr(i, kmerSize);
        filter.insert(temp);
        size_B_est++;
    }

    int int_est = 0;
    for (int i = 0; i <= strings[1].size() - kmerSize; i++) {
        temp = strings[1].substr(i, kmerSize);
        if(filter.contains(temp)) {
            int_est++;
        }
        size_A++;
    }

    int_est -= floor(FALSE_POSTIVITY * 10);  // adjust for the false positive rate
    float containment_est = int_est / float(10);  //estimate of the containment index
    float jaccard_est = (size_A * containment_est) / (float)(size_A + size_B_est - (size_A * containment_est));
    cout << "\nJaccard similarity between index using containment hash approach is: "<<jaccard_est<<"\n";

}