#include "HepLib.h"
#include <dirent.h>

using namespace HepLib;

int main(int argc, char** argv) {
        
    DIR *di;
    struct dirent *dir;
    di = opendir(".");
    set<string> pns;
    if (di) {
        while ((dir = readdir(di)) != NULL) {
            string fn = dir->d_name;
            if(string_end_with(fn, ".start")) {
                string_replace_all(fn, ".start", ".tables");
                if(!file_exists(fn)) pns.insert(fn);
            }
        }
        closedir(di);
    }
    if(pns.size()>0) {
        cout << pns.size() << " tables NOT exist: " << endl;   
        for(auto pn : pns) cout << pn << " ";
        cout << endl;
    } else cout << "ALL tables exist." << endl;
    return 0;
}
