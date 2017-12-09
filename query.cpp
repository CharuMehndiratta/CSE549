#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include "MurmurHash3.h"

using namespace std;


void generate_sketch() {
    
}

void genertae_jacard_index(string long_read) {
    

}

void read_dataset(string filename) {

    string sequence, line;
    ifstream file (filename);

    if (file.is_open()) {
        getline(file, line);
        while (getline(file, line)) {
            if (line[0] != '>') {

                genertae_jacard_index(line);

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
    char *long_read_file = NULL, *bloom_filter_file = NULL, *min_sketch_file = NULL;

    /* Get options from command line arguments */
    while ((option = getopt(argc, argv, "r:b:m")) != -1) {
        switch (option) {
            case 'r':
                long_read_file = optarg;
                break;
            case 'b':
                bloom_filter_file = optarg;
                break;
            case 'm':
                min_sketch_file = optarg;
                break;
            default:
                exit(EXIT_FAILURE);
        }
    }

    if (long_read_file == NULL) {

    }

    read_dataset(long_read_file);

    return 0;
}
