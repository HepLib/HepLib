#include "HepLib.h"
#include <dirent.h>

using namespace HepLib;

int main(int argc, char** argv) {
    
    if(argc<3) {
        cout << "usage: " << argv[0] << " <ostr> <nstr>" << endl;
        return 0;
    }
    
    string o_str = argv[1];
    string_replace_all(o_str, "\\n", "\n");
    string n_str = argv[2];
    string_replace_all(n_str, "\\n", "\n");
    DIR *di;
    struct dirent *dir;
    di = opendir(".");
    if (di) {
        while ((dir = readdir(di)) != NULL) {
            string fn = dir->d_name;
            if(string_end_with(fn, ".config")) {
                cout << "\r                                             \r" << flush;
                cout << " - Replacing " << fn << " ..." << flush;
                string ostr = file2str(fn);
                if(n_str.find("<n>")==string::npos) string_replace_all(ostr, o_str, n_str);
                else {
                    string fn2 = fn.substr(0, fn.find_last_of("."));
                    string n_str2 = n_str;
                    string_replace_all(n_str2, "<n>", fn2);
                    string_replace_all(ostr, o_str, n_str2);
                }
                str2file(ostr, fn);
            }
        }
        closedir(di);
    }    
    cout << endl;
    return 0;
}
