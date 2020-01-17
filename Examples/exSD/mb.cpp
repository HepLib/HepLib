#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    auto vs = SD::vs;
    
    if(false) {
        XIntegrand xint;
        xint.Functions = lst{ 1, vs*x(0)-3*x(2) };
        xint.Exponents = lst{ 1, -1+ep };
            
        SD work;
        work.epN = 0;
        work.Verbose = 12;
        work.sN = 300;
        
        if(work.SecDec==NULL) work.SecDec = new SecDecG();
        if(work.Integrator==NULL) work.Integrator = new HCubature();
        if(work.Minimizer==NULL) work.Minimizer = new MinUit();
        
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "I*(1.0471975511965843 +- 2.94385E-7)+(-0.6365141682948453 +- 9.6122E-8)" << endl;
    }
    
    if(true) {
        XIntegrand xint;
        ex r = ex(1)/3;
        xint.Functions = lst{ 1, x(0)-vs*x(2) };
        xint.Exponents = lst{ 1, -1+ep };
            
        SD work;
        work.epN = 0;
        work.Verbose = 12;
        work.sN = 30;
        work.CheckEnd = true;
        
        if(work.SecDec==NULL) work.SecDec = new SecDecG();
        if(work.Integrator==NULL) work.Integrator = new HCubature();
        if(work.Minimizer==NULL) work.Minimizer = new MinUit();
        
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
        cout << "check with:" << endl;
        cout << "I*(1.0472006821272757 +- 1.66998E-5)+(-0.26162451065437947 +- 1.70888E-5)" << endl;
    }
    
    return 0;
}
