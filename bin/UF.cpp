#include "SD.h"

using namespace HepLib;

auto ep = SD::ep;
auto eps = SD::eps;
auto iEpsilon = SD::iEpsilon;

symtab table = { {"ep", ep}, {"eps", eps}, {"iEpsilon", iEpsilon} };
parser reader(table);

ex exRead() {
    ostringstream oss;
    while(true) {
        string line;
        getline(cin, line);
        if(line == ".end") break;
        oss << line << endl;
    }
    ex res = reader(oss.str().c_str());
    return res;
}

int main(int argc, char** argv) {
        
    ex ls = exRead();
    ex tls = exRead();
    ex ps = exRead();
    ex ns = exRead();
    
    ex lr = exRead();
    ex tlr = exRead();
    ex nr = exRead();
    
    FeynmanParameter fp;
    
    fp.LoopMomenta = ex_to<lst>(ls);
    fp.tLoopMomenta = ex_to<lst>(tls);
    fp.Propagators = ex_to<lst>(ps);
    fp.Exponents = ex_to<lst>(ns);
    
    for(auto vv : lr) {
        fp.lReplacements[vv.op(0)] = vv.op(1);
    }
    
    for(auto vv : tlr) {
        fp.tReplacements[vv.op(0)] = vv.op(1);
    }
    
    for(auto vv : nr) {
        fp.nReplacements[vv.op(0)] = vv.op(1);
    }
    
    fp.nReplacements[ep] = ex(1)/11;
    fp.nReplacements[eps] = ex(1)/111;
    
    SD work;
    work.Verbose = 0;
    //work.CheckF1 = true;
    
    work.Initialize(fp);
    
    cout << "{";
    int tot = work.FunExp.size();
    int cur = 1;
    for(auto kv : work.FunExp) {
        cout << "{" << kv.first << ", " << kv.second << "}";
        if(cur<tot) cout << ", ";
        cur++;
    }
    cout << "}";
            
    return 0;
}