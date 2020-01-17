#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {

    //SD::debug = true;
    
    auto ep = SD::ep;
    auto eps = SD::eps;
    auto iep = SD::iEpsilon;
    auto t = SD::vs;
    
    Digits = 30;
        
    XIntegrand xint;
    
    xint.Functions = {x(1), t*x(2)};
    xint.Exponents = {1+ep, 1+ep};
    
    SD work;
    work.epN = 1;
    work.sN = 4;
    work.Verbose = 12;
    //work.ParallelProcess = 0;
    
    //work.use_cpp = false;
    //SD::debug = true;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    work.Evaluate(xint);

    work.VEPrint();
    
    cout << WHITE << "Check with: " << RESET << endl;
    cout << "(s*(0.250000000000000000004079635449552256 +- 4.6192E-21)) + (s*(-0.2499999004825011 +- 1.2666E-5)+log(s)*s*(0.250000000000000000004079635449552256 +- 4.6192E-21))*ep" << endl;
        
    return 0;
}
