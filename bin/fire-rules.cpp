#include "HepLib.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    if(argc!=2 && argc!=3) {
        cout << "usage: " << argv[0] << " <problem number> [o], o for old version compatable." << endl;
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
    
    if(argc==3) { // to support old versions
        lst mul_repl;
        for(auto const & item : fire.Replacement) {
            if(is_a<mul>(item.op(0))) {
                bool found = false;
                for(auto const & it : item) if(is_a<wildcard>(it)) { found = true; break; }
                if(!found) mul_repl.append(wild(10000)*item.op(0) == wild(10000)*item.op(1));
            }
        }
        for(auto const & item : mul_repl) fire.Replacement.append(item);
    }
    fire.Integral.append(0);
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
