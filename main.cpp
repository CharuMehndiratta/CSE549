#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <unordered_set>

#define MAX_PRIME 9999999999971
#define MAX_HASH_ALGO 5

using namespace std;

float jaccard_similarity(vector<int> & read1, vector<int>& read2){
    unordered_set<int> set1, set2;
    float common = 0, total = 0;
    float jaccard;
    
    for(int i = 0; i < read1.size(); i++){
        set1.insert(read1[i]);
    }
    
    for(int j = 0; j < read2.size(); j++){
        set2.insert(read2[j]);
    }
    total = set1.size();
    
    for(auto it = set2.begin(); it != set2.end(); it++){
        if(set1.find(*it) != set1.end())
            common++;
        else
            total++;
    }
    cout << common << " " << total <<"\n";
    jaccard = common/total;
    return jaccard;
}


// https://stackoverflow.com/questions/15518418/whats-behind-the-hashcode-method-for-string-in-java
int random_hash(string str) {
    int hash = 0, multiplier = 1, count = int(str.size());
    for (int i = count - 1; i >= 0; i--) {
        hash += str[i] * multiplier;
        int shifted = multiplier << 5;
        multiplier = shifted - multiplier;
    }
    return hash;
}

vector<int> random_hash_vector;
void generate_rest_random_hash() {
    srand (time(NULL));
    for(int i=0;i<MAX_HASH_ALGO;i++) {
        cout<<"\n rand -- "<<rand();
        random_hash_vector.push_back(rand());
    }
    
}

void generate_hash_value(int k, string str, vector<int> &hash_value, vector<string> &hash_string) {
    
    int base_hash, temp_hash;
    string temp;
    cout<<"\n input  "<<str;
    cout<<"\n";
    for (int i = 0; i <= str.size()-k; i++) {
        temp = str.substr(i,k);
        cout<<"\n temp  "<<temp;
        base_hash = random_hash(temp) % MAX_PRIME;
        if(hash_value[0] > base_hash) {
            hash_value[0] = base_hash;
            hash_string[0] = temp;
        }
        cout<<" hash -- "<<base_hash;
        for( int j = 1; j < MAX_HASH_ALGO; j++) {
            temp_hash = ( random_hash_vector[j] ^ base_hash) % MAX_PRIME;
            cout<<" hash -- "<<temp_hash;
            if(hash_value[j] > temp_hash) {
                hash_value[j] = temp_hash;
                hash_string[j] = temp;
            }
        }
    }
}


int main() {
    string input = "", tmp;
    vector < vector<int> > hash_values;
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
        vector<int> tmp(200, INT_MAX);
        vector<string> temp_str(200);
        hash_values.push_back(tmp);
        hash_strings.push_back(temp_str);
        generate_hash_value(3, strings[i], hash_values[i], hash_strings[i]);
        
        cout<<"-------------\n";
        for(int j = 0; j<MAX_HASH_ALGO; j++) {
            cout<<hash_values[i][j]<<" => ";
            cout<<hash_strings[i][j]<<" ";
        }
    }
    
    for(int k = 0; k < hash_values.size()-1; k++){
        for(int l = k+1; l < hash_values.size(); l++){
            cout << "Jaccard similarity between " << k << " and " << l << " is: ";
            cout << jaccard_similarity(hash_values[k], hash_values[l]);
            cout << "\n";
        }
    }
    
}

