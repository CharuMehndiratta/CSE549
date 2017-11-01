#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <unordered_set>
#include "MurmurHash3.h"
#include <string>

#define MAX_HASH_ALGO 5
#define MAX_HASH_SIZE 10000000000

using namespace std;

float jaccard_similarity(vector<uint64_t> & read1, vector<uint64_t>& read2){
    unordered_set<uint64_t> set1, set2;
    float common = 0, total = 0;
    float jaccard;
    
    for(int i = 0; i < read1.size(); i++){
        set1.insert(read1[i]);
    }
    
    for(int j = 0; j < read2.size(); j++){
        set2.insert(read2[j]);
    }
    total = set1.size();

    unordered_set<uint64_t> ::iterator it;
    
    for(it = set2.begin(); it != set2.end(); it++){
        if(set1.find(*it) != set1.end())
            common++;
        else
            total++;
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
        cout<<"\n rand -- "<<x;
        random_hash_vector.push_back(x);
    }

}


void generate_hash_value(int k, string str, vector<uint64_t> &hash_value, vector<string> &hash_string) {
    
    uint64_t hash_otpt[2], temp_hash;
    string temp;
    cout<<"\n input  "<<str;
    cout<<"\n";
    for (int i = 0; i <= str.size()-k; i++) {
        temp = str.substr(i,k);
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
    string input = "", tmp;
    vector < vector<uint64_t> > hash_values;
    vector < vector<string> > hash_strings;
    vector < string > strings;
    int n;
    ifstream file ("dataset.txt");
    
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
        vector<uint64_t> tmp(200, LLONG_MAX);
        vector<string> temp_str(200);
        hash_values.push_back(tmp);
        hash_strings.push_back(temp_str);
        generate_hash_value(3, strings[i], hash_values[i], hash_strings[i]);
        
        cout<<"\n-------------\n";
        for(int j = 0; j<MAX_HASH_ALGO; j++) {
            cout<<hash_values[i][j]<<" => ";
            cout<<hash_strings[i][j]<<" ";
        }
    }
    
    for(int k = 0; k < hash_values.size()-1; k++){
        for(int l = k+1; l < hash_values.size(); l++){
            cout << "\nJaccard similarity between index " << k << " and index " << l << " is: ";
            cout << jaccard_similarity(hash_values[k], hash_values[l]);
            cout << "\n";
        }
    }
    
}