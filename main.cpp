#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>

#define MAX_PRIME 9999999999971
#define MAX_HASH_ALGO 5

using namespace std;

// https://stackoverflow.com/questions/15518418/whats-behind-the-hashcode-method-for-string-in-java
int random_hash(string str) {
    int hash = 0, multiplier = 1, count = str.size();
    for (int i = count - 1; i >= 0; i--) {
        hash += value[i] * multiplier;
        int shifted = multiplier << 5;
        multiplier = shifted - multiplier;
    }
    return hash;
}

void generate_hash_value(int k, string str, vector < vector<int> > &hash_value, vector < vector<int> > &hash_string) {
	srand (time(NULL));
	int base_hash, temp_hash;
	string temp;
	int min = INT_MAX;

	for (int i = 0; i < str.size(); i = i + k) {
		temp = str.substr(0,k);
		base_hash = random_hash(temp) % MAX_PRIME;
		for( int j = 0; j < MAX_HASH_ALGO; j++) {
			temp_hash = (rand() ^ base_hash) % MAX_PRIME;
			if(hash_value[j] > temp_hash) {
				hash_value[j] = temp_hash;
				hash_string[j] = temp;
			}
		}
	}
}


// void hash_shingles(string first, string second, vector<int> first_hash, vector<string> first_string_hash, vector<int> second_hash, vector<string> second_string_hash) {
// 	string temp;

// 	for (int i = 0; i < first.size(); i = i + 3) {
// 		temp = str.substr(0, 3);
// 		if (i == 0) {
// 			generate_hash_value(temp);
// 		}
		
		
// 	}
// }

int main() {

	string first, second;
	// ifstream file ("rosalind_dataset.txt");

    // if (file.is_open()) {
    // 	getline(file, input);
    // 	while(getline(file, input) && input[0] != '>') {
    //         first += input;
    //     }
    //     while(getline(file, input)) {
    //         second += input;
    //     }
    // } else {
    //     cout<<"\n Unable to open file";
    // }
    // file.close();

    first = "ACGTACGTCCAT";
    second = "ACAAGCTACTCGTACGGCCAT";
    int n = 2;

    vector < vector<int> > hash_values;
    vector < vector<string> > hash_strings;
    vector < vector<string> > strings;

    // for (int i = 0; i < n; i++) {
    // 	strings[i].push_back()
    // }

    strings[0].push_back(first);
    strings[1].push_back(second);


    for(int i = 0; i < n; i++) {
    	vector<int> tmp(200, INT_MAX);
    	hash_values.push_back(tmp);
    	generate_hash_value(strings[i], hash_values[i], hash_strings[i]);
    }
    










}