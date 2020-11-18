#include "BASIC.h"
#include <getopt.h>

using namespace std;
using namespace HepLib;

int main(int argc, char** argv) {

    if(argc<2) {
        cout << "Usage: " << argv[0] << " [--prefix|--cflags|--ldflags|--help]" << endl;
        return 0;
    }
        
    // handle long options
    int opt;
    int digit_optind = 0;
    int option_index = 0;
    const char *string = "";
    static struct option long_options[] = {
        { "prefix", no_argument, NULL, 'p' },
        { "cflags", no_argument, NULL, 'c' },
        { "ldflags", no_argument, NULL, 'l' },
        { "help", no_argument, NULL, 'h' },
        {NULL, 0, NULL, 0},
    };
    while((opt =getopt_long_only(argc,argv,string,long_options,&option_index))!= -1) {
        switch (opt) {
            case 'p': cout << InstallPrefix; break;
            case 'c': cout << INC_FLAGS; break;
            case 'l': cout << LIB_FLAGS; break;
            default: cout << "Usage: " << argv[0] << " [--prefix|--cflags|--ldflags|--help]" << endl; break;
        }
    }

    return 0;
}
