#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {

    //SD::debug = true;
    
    auto ep = SD::ep;
    auto eps = SD::eps;
    auto iep = SD::iEpsilon;
    auto t = SD::vs;
            
    XIntegrand xint;
    
    xint.Functions = lst{x(1), x(1)+t*x(2)};
    xint.Exponents = lst{1+2*ep, 1+ep};
    
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
    cout << "((0.3333333333 +- 9.369670984E-48)+(0 +- 5.041104112E-49)*s^3+(0.25 +- 3.009265538E-35)*s) + ((-0.33333333304091733 +- 6.1011E-6)+(-0.04513888888803162 +- 1.49863E-9)*s^3+(0.1666666667 +- 4.684835492E-48)*s^2+s^4*(-0.01666666667 +- 1.54751922E-49)+(0.04166666667 +- 1.069379671E-48)*log(s)*s^3+s*(-0.1249999302627286 +- 8.7726E-6))*ep" << endl;
        
    return 0;
}
