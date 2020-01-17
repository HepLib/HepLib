#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    symbol q1("q1"), q2("q2");
    symbol p1("p1"), k1("k1");
    symbol m("m"), s("s");
    
    ex m2=m*m;
    FeynmanParameter fp;
    
    fp.LoopMomenta = lst{q1,q2};
    fp.Propagators = lst {
-pow(-p1 + q1,2), -pow(p1 + q1,2), m2-pow(q2,2), m2-pow(-k1 + q2,2), m2-pow(-p1 + q1 + q2,2), m2-pow(-k1 + p1 + q1 + q2,2)
    };

    fp.Exponents = lst{ 1, 1, 1, 1, 1, 1 };
    fp.lReplacements[p1*p1] = m2;
    fp.lReplacements[k1*k1] = 0;
    fp.lReplacements[k1*p1] = m2;
    fp.lReplacements[m] = 1;
    fp.Prefactor = pow(I*pow(Pi,2-ep)*exp(-ep*Euler), -2);
    
    SD work;
    work.epN = 0;
    work.Verbose = 100;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    work.Evaluate(fp);

    cout << work.VEResult() << endl;
    cout << "check with:" << endl;
    cout << "(-0.027330019 +- 4.30844E-5)+I*(0.37515628 +- 4.3997E-5)" << endl;
    
    return 0;
}
