#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include "MurmurHash3.h"
// #include "BloomFilter.hpp"

using namespace std;

/* Number of hash functions, false positive and kmer size with default values */
double false_positive = 0.001;
int kmer_size = 16;
int num_hash = 100;

/*************************************************************/
/* Start execution                                           */
/*************************************************************/
int main(int argc, char *argv[]) {
    int option;
    char *reference_file = NULL;

    /* Get options from command line arguments */
    while ((option = getopt(argc, argv, "r:f:k:h:")) != -1) {
        switch (option) {
            case 'r':
                reference_file = optarg;
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

    return 0;
}
