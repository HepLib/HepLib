#include "HepLib.h"

using namespace HepLib;
using namespace SD;

int main(int argc, char** argv) {

    Symbol k("k"), l("l"), r("r"), q("q"), q1("q1"), q2("q2");
    Symbol p1("p1"), p2("p2"), p3("p3"), p4("p4"), p5("p5");
    Symbol m("m"), s("s");
    #define t vs

    FeynmanParameter fp;
    
    fp.LoopMomenta = lst {k};
    fp.Propagator = lst{ -pow(k,2),-pow(k + p1,2),-pow(k + p1 + p2,2),-pow(k + p1 + p2 + p4,2) };
    fp.Exponent = lst{ 1, 1, 1, 1 };
    fp.lReplacement[p1*p1] = 0;
    fp.lReplacement[p2*p2] = 0;
    fp.lReplacement[p4*p4] = 0;
    fp.lReplacement[p1*p2] = -s/2;
    fp.lReplacement[p2*p4] = -t/2;
    fp.lReplacement[p1*p4] = (s+t)/2;
    fp.lReplacement[s] = 1;
    fp.Prefactor = pow(I*pow(Pi,2-ep)*exp(-ep*Euler), -1);
        
    SecDec work;
    work.eps_lst = lst{ lst{ep, 1} };
    work.vsN = 1;
    Verbose = 100;
    
    auto intor = new HCubature();
    
    intor->QXLimit = 1E-1;
    intor->MPXLimit = 1E-3;
    
    work.Integrator = intor;
            
    work.Evaluate(fp);

    cout << work.VEResult() << endl;
    
    return 0;
}
