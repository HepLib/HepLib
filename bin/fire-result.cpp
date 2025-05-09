#include "HepLib.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    if(argc!=2 && argc != 3) {
        cout << "usage: " << argv[0] << " <problem number> [1: MI|2: IBP rules|3: MI rules]" << endl;
        return 0;
    }
    
    string pn = argv[1];
    string opt = "";
    if(argc==3) opt = argv[2];
    
    string garfn = pn+".gar";
    if(!file_exists(garfn)) {
        cout << "The file: " << garfn << " NOT found." << endl;
        return 0;
    }
    FIRE fire;
    fire.FROM(garRead(garfn));
    fire.WorkingDir = ".";
    fire.Integral.append(0);
    fire.Import();
    auto rm = fire.FindRules(true);
    
    lst mi_rules;
    for(auto const & kv : rm.first) {
        mi_rules.append(kv.first == kv.second);
    }
    
    if(argc==3) {
        if(opt=="1") cout << fire.MIntegral << endl;
        else if(opt=="2") cout << fire.Rules << endl;
        else if(opt=="3") cout << mi_rules << endl;
        else cout << "opt: " << opt << " Not supported." << endl;
    } else {
        cout << lst{ fire.MIntegral, fire.Rules, mi_rules } << endl;
    }
    
    return 0;
}
