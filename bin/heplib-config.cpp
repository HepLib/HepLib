#include "BASIC.h"
#include <getopt.h>

using namespace std;
using namespace HepLib;

int main(int argc, char** argv) {
        
    // handle long options
    int opt;
    int digit_optind = 0;
    int option_index = 0;
    const char *string = "";
    static struct option long_options[] = {
        { "prefix", no_argument, NULL, 'p' },
        { "cflag", no_argument, NULL, 'c' },
        { "ldflag", no_argument, NULL, 'l' },
        { "help", no_argument, NULL, 'h' },
        {NULL, 0, NULL, 0},
    };
    while((opt =getopt_long_only(argc,argv,string,long_options,&option_index))!= -1) {
        switch (opt) {
            case 'p': cout << InstallPrefix; break;
            case 'c': cout << INC_FLAGS; break;
            case 'l': cout << LIB_FLAGS; break;
            default: cout << "heplib-config [--prefix] [--cflag] [--ldflag] [--help]" << endl; break;
        }
    }
    

    return 0;
}
