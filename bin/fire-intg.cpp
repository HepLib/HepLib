#include "HepLib.h"
#include <dirent.h>

using namespace HepLib;

int main(int argc, char** argv) {
    
    if(argc<3) {
        cout << "usage: " << argv[0] << " -s / -c pn" << endl;
        cout << "       -s to split pn.intg to pn-i.intg, including pn.config" << endl;
        cout << "       -c to collect pn-i.intg to pn.intg" << endl;
        return 0;
    }
    
    string act = argv[1];
    string pn = argv[2];
    
    if(act=="-s") {
        string fn_intg = pn+".intg";
        if(!file_exists(fn_intg)) {
            cout << "the integrals file: " << fn_intg << " NOT found." << endl;
            return 0;
        }
        auto intg = file2ex(fn_intg);
    } else if(act=="-c") {
    
    
    }
    return 0;
}
