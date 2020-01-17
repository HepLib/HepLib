#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    {
        XIntegrand xint;
        xint.Functions = lst{ x(1), x(2) };
        xint.Exponents = lst{ -1+ep, -1+ep };
        xint.Deltas.push_back(lst{x(1),x(2)});
        
        SD work;
        work.epN = 3;
        work.Verbose = 100;
        
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        
        work.Evaluate(xint);
        
        cout << work.VEResult() << endl;
    }
    
    {
        XIntegrand xint;
        xint.Functions = lst{ x(1), x(2), x(1), x(0)+x(1), x(1)+x(2) };
        xint.Exponents = lst{ -1+ep, -1+ep, 1, -2, -2*ep };
        xint.Deltas.push_back(lst{x(0), x(1),x(2)});
        
        SD work;
        work.epN = 3;
        work.Verbose = 100;
        
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        
        work.Evaluate(xint);
        
        cout << work.VEResult() << endl;
    
    }
    
    return 0;
}
