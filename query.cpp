#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include "MurmurHash3.h"

using namespace std;

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

    return 0;
}
