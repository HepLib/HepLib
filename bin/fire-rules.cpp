#include "HepLib.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    if(argc!=2) {
        cout << "usage: " << argv[0] << " <problem number>" << endl;
        return 0;
    }
    
    string pn = argv[1];
    
    string garfn = pn+".gar";
    if(!file_exists(garfn)) {
        cout << "The file: " << garfn << " NOT found." << endl;
        return 0;
    }
    FIRE fire;
    fire.FROM(garRead(garfn));
    fire.WorkingDir = ".";
    fire.Integrals.append(0);
    fire.Import();
    auto rm = fire.FindRules(true);
    ofstream rules_out("./"+pn+".rules");
    for(auto const & r : rm.first) {
        ostringstream oss;
        oss << ex2str(r.first.subs(F(w1,w2)==F(0,w2)));
        oss << " -> ";
        oss << "{{1," << ex2str(r.second.subs(F(w1,w2)==F(0,w2))) << "}}";
        string s = oss.str();
        string_replace_all(s, "F(", "G(");
        string_replace_all(s, ",", ", ");
        rules_out << s << endl << endl;
    }
    rules_out.close();
    return 0;
}
