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
        map<ex,lst,ex_is_less> sector_intgs;
        for(const auto & i : intg) {
            lst key;
            for(const auto & ii : i.op(1)) {
                if(ii>0) key.append(1);
                //else if(ii<0) key.append(-1);
                else key.append(0);
            }
            sector_intgs[key].append(i);
        }
        int n = 0;
        string config = file2str(pn+".config");
        for(auto const & kv : sector_intgs) {
            n++;
            string ns = to_string(n);
            string cc = config;
            string_replace_all(cc, "#database db"+pn, "#database db"+pn+"-"+ns);
            string_replace_all(cc, "#integrals "+pn+".intg", "#integrals "+pn+"-"+ns+".intg");
            string_replace_all(cc, "#output "+pn+".tables", "#output "+pn+"-"+ns+".tables");
            str2file(cc, pn+"-"+ns+".config");
            ex2file(kv.second, pn+"-"+ns+".intg");
        }
    } else if(act=="-c") {
    
    
    }
    return 0;
}

