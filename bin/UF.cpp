#include "SD.h"

using namespace HepLib::SD;

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
    
    ex isQuasi = exRead();
    ex ls = exRead();
    ex tls = exRead();
    ex ps = exRead();
    ex ns = exRead();
    
    ex lr = exRead();
    ex tlr = exRead();
    ex nr = exRead();
    
    FeynmanParameter fp;
    
    if(isQuasi>0) fp.isQuasi = true;
    fp.LoopMomenta = ex_to<lst>(ls);
    fp.tLoopMomenta = ex_to<lst>(tls);
    fp.Propagator = ex_to<lst>(ps);
    fp.Exponent = ex_to<lst>(ns);
    
    for(auto vv : lr) {
        fp.lReplacement[vv.op(0)] = vv.op(1);
    }
    
    for(auto vv : tlr) {
        fp.tReplacement[vv.op(0)] = vv.op(1);
    }
    
    for(auto vv : nr) {
        fp.nReplacement[vv.op(0)] = vv.op(1);
    }
    
    fp.nReplacement[ep] = ex(1)/11;
    fp.nReplacement[eps] = ex(1)/111;
    
    SecDec work;
    //work.CheckEnd = true;
    
    work.Initialize(fp);
    
    cout << "{";
    int tot = work.FunExp.size();
    int cur = 1;
    for(auto fe : work.FunExp) {
        cout << "{" << fe.op(0) << ", " << fe.op(1) << "}";
        if(cur<tot) cout << ", ";
        cur++;
    }
    cout << "}";
            
    return 0;
}
