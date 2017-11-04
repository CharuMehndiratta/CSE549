#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <string>
#include <cmath>
#include "MurmurHash3.h"
#include "BloomFilter.hpp"
#include <unordered_set>
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

   // cout << common << " / " << total <<"\n";
    jaccard = common/total;
    cout << jaccard << "  ";
    return jaccard;
}


vector<uint64_t> random_hash_vector;
void generate_rest_random_hash() {
    srand (time(NULL));
    for (int i = 0; i < MAX_HASH_ALGO; i++) {
        int x = rand();
       // cout<<"\n rand -- "<<i;
        random_hash_vector.push_back(i);
    }

}


void generate_hash_value(int k, string str, vector<uint64_t> &hash_value, vector<string> &hash_string) {
    
    uint64_t hash_otpt[2], temp_hash;
    string temp;
 //   cout<<"\n input  "<<str;
   // cout<<"\n";
    for (int i = 0; i <= str.size() - k; i++) {
        temp = str.substr(i, k);
        const char *key = temp.c_str();
     //   cout<<"\n key - "<<temp;
        for( int j = 0; j < MAX_HASH_ALGO; j++) {
            MurmurHash3_x64_128(key, (uint64_t)strlen(key), random_hash_vector[j], hash_otpt);
            temp_hash = (hash_otpt[0] +  MAX_HASH_ALGO * hash_otpt[1]) % MAX_HASH_SIZE;
       //     cout<<" temp_hash "<<temp_hash;

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
    vector < vector<uint64_t> > hash_reads;
    vector < vector<uint64_t> > hash_references;
    vector < vector<string> > hash_kmer_reads;
    vector < vector<string> > hash_kmer_references;
    vector < string > long_reads;
    vector<string> references;
    vector<float> min_hash_jaccard;
    vector<float> con_hash_jaccard;
    vector<float> true_jaccard;

    int n;
    int kmerSize;
    cout << "Enter KmerSize: ";
    cin >> kmerSize;
    cout << "\n";

    /*************************************************************/
    /* Read the long_reads from file	 			             */
    /*************************************************************/
    ifstream file ("long_reads.txt");
    if (file.is_open()) {
        getline(file, tmp);
        while(getline(file, tmp)) {
            if(tmp[0] != '>') {
                input += tmp;
            }
            else {
                long_reads.push_back(input);
                input = "";
            }
        }
	long_reads.push_back(input);
    } else {
        cout<<"\n Unable to open Long_Reads file";
    }
    file.close();
    
    /*************************************************************/
    /* Read the references from file	 			             */
    /*************************************************************/
    tmp = "";
    input = "";
    ifstream file2 ("references.txt");
    if (file2.is_open()) {
	 getline(file2, tmp);
	 while(getline(file2, tmp)) {
             if(tmp[0] != '>') {
                  input += tmp;
	     }
             else {
                   references.push_back(input);
                   input = "";
             }
	 }
        references.push_back(input);
    } 
    else {
        cout<<"\n Unable to open References file";
    }
    file2.close();

    generate_rest_random_hash();
    /*************************************************************/
    /* Initialize the min_hash values and kmer for long reads 	 */
    /*************************************************************/
    n = int(long_reads.size());
    for(int i = 0; i < n; i++) {
        vector<uint64_t> tmp(MAX_HASH_ALGO, LLONG_MAX);
        vector<string> temp_str(MAX_HASH_ALGO);
        hash_reads.push_back(tmp);
        hash_kmer_reads.push_back(temp_str);
        generate_hash_value(kmerSize, long_reads[i], hash_reads[i], hash_kmer_reads[i]);
    }
    n = 0;

    /*************************************************************/
    /* Initialize the min_hash values and kmer for references 	 */
    /*************************************************************/
    n = int(references.size());
    for(int i = 0; i < n; i++) {
        vector<uint64_t> tmp(MAX_HASH_ALGO, LLONG_MAX);
        vector<string> temp_str(MAX_HASH_ALGO);
        hash_references.push_back(tmp);
        hash_kmer_references.push_back(temp_str);
        generate_hash_value(kmerSize, references[i], hash_references[i], hash_kmer_references[i]);
    }
    
    /*************************************************************/
    /* Generate Jaccard similarity between long reads		     */
    /* and references using min_hash			 	             */
    /*************************************************************/
    for(int k = 0; k < hash_references.size(); k++){
        for(int l = 0; l < hash_reads.size(); l++){
	    min_hash_jaccard.push_back(jaccard_similarity(hash_references[k], hash_reads[l]));
        //cout << min_hash_jaccard.back();
        }
    }

    /*************************************************************/
    /* Generate Jaccard similarity between long reads		     */
    /* and references using containment_hash		 	         */
    /*************************************************************/
    for(int i = 0; i < references.size(); i++){
	
        int capacity = 2 * references[i].size();
        bloom_parameters parameters;
    
        // How many elements roughly do we expect to insert?
        parameters.projected_element_count = capacity;

        // Maximum tolerable false positive probability? (0,1)
        parameters.false_positive_probability = FALSE_POSTIVITY;

        // Simple randomizer (optional)
        parameters.random_seed = rand();

        //cout<<"\n parameters.random_seed  "<<parameters.random_seed ;
        if (!parameters) {
            cout << "Error - Invalid set of bloom filter parameters!";
            return 1;
        }

        parameters.compute_optimal_parameters();
        bloom_filter filter(parameters);

        string temp;
        int size_B_est = 0; // String 0
        int size_A = 0; // String 1

        //Inserting references[i] in Bloom Filter
        for (int j = 0; j <= references[i].size() - kmerSize; j++) {
            temp = references[i].substr(j, kmerSize);
            if(!filter.contains(temp)){
                filter.insert(temp);
                size_B_est++;
            }
        }
        
        int int_est = 0;
        int hash_value;
        for(int k = 0; k < long_reads.size(); k++){
	 	
            unordered_set<string> read_set;
         	for (int l = 0; l <= long_reads[k].size() - kmerSize; l++) {
             	  temp = long_reads[k].substr(l, kmerSize);
                  if(filter.contains(temp)) {
                      int_est++;
                   }
                   read_set.insert(temp);
                }

            hash_value = references[i].size()/long_reads[k].size();
            size_A = read_set.size();
         	int_est -= floor(FALSE_POSTIVITY * hash_value);  // adjust for the false positive rate
	        float containment_est = int_est / float(hash_value);  //estimate of the containment index
       		float jaccard_est = (size_A * containment_est) / (float)(size_A + size_B_est - (size_A * containment_est));
            con_hash_jaccard.push_back(jaccard_est);
            //cout <<  con_hash_jaccard.back();
      	 }   	
     }

    /*************************************************************/
    /* Generate true Jaccard similarity between long reads	     */
    /* and references using conventional method		 	         */
    /*************************************************************/

    for(int i = 0; i < references.size(); i++){
        unordered_set<string> ref_kmers;
        for(int j = 0; j < references[i].size() - kmerSize; j++){
            string  temp = references[i].substr(j, kmerSize);
            ref_kmers.insert(temp);
        }
        for(int k = 0; k < long_reads.size(); k++){
            unordered_set<string> read_kmers;
            for(int l = 0; l < long_reads[k].size() - kmerSize; l++){
                string tmp =  long_reads[k].substr(l, kmerSize);
                read_kmers.insert(tmp);
            }
            
            float un = ref_kmers.size();
            float intr = 0;
            
            for(auto it = read_kmers.begin(); it != read_kmers.end(); it++){
                if(ref_kmers.find(*it) != ref_kmers.end())
                    intr++;
                else
                    un++;
            }
            true_jaccard.push_back(intr/un);
        }
    }

    cout << "\n";
    cout << "Reference Idx \t" << "Long Read Idx \t" << "Min_Hash  \t" << "Containment_Hash \t" << "True Jaccard\n";
    int len = 0;
    for(int i = 0; i < references.size(); i++){
        for(int j = 0; j < long_reads.size(); j++){
            cout << i << "\t\t" << j << "\t\t" << min_hash_jaccard[len] << "\t\t" << con_hash_jaccard[len] << "\t\t" << true_jaccard[len] << "\n";
            len++;
        }
    }
}	
