#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    symbol m("m");
    
    ex m2=m*m;
    XIntegrand xint;
    
    xint.Functions = lst{ x(4)*(x(1) + 2*x(4) + 4*x(4)*x(5) + 2*x(1)*x(6) + 4*x(4)*x(6)),
     (x(4)*(x(1) + 2*x(4) + 4*x(4)*x(5) + 2*x(1)*x(6) + 4*x(4)*x(6))*
        (3*x(4) + x(2)*x(5) + 6*x(4)*x(5) + 3*x(1)*x(6) + 6*x(4)*x(6))),32*x(5)*x(6),48*x(5)*x(6),1,
     x(5)*x(6) };
    xint.Exponents = lst{ 3*ep,-1 - 2*ep,-3*ep,1 + 2*ep,1,-1 };
    xint.nReplacements[ep] = ex(1)/11;
    
    xint.Deltas.push_back(lst{x(1),x(2),x(4)});
    xint.Deltas.push_back(lst{x(5),x(6)});
    
    SD work;
    work.epN = -1;
    work.Verbose = 2;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    work.Evaluate(xint);
    
    cout << work.VEResult() << endl;
    
    return 0;
}
