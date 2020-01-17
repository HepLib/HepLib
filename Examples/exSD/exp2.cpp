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
    
    xint.Functions = lst{x(1), x(1)+t*x(2), 1+t, 1-2*t};
    xint.Exponents = lst{0, -3+ep, -3, 3};
    
    SD work;
    work.epN = 2;
    work.sN = 3;
    work.Verbose = 12;
    work.ParallelProcess = 1;
    //work.ParallelProcess = 0;
    
    //work.use_cpp = false;
    //SD::debug = true;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    work.Evaluate(xint);

    work.VEPrint();
    
    cout << WHITE << "Check with: " << RESET << endl;
    cout << "((207.5 +- 1.245493417E-45)*s^3+(50 +- 3.414424874E-47)*s+(-18.5 +- 0)+s^(-2)*(-0.5 +- 0)+(4.5 +- 0)*s^(-1)+(-108.5 +- 3.081022244E-46)*s^2) + (log(s)*(45 +- 0)*s+(-0.75 +- 0)*s^(-2)+log(s)*s^3*(139.5 +- 0)+(-27.25 +- 0)+s^3*(229.9583333 +- 1.747670738E-45)+log(s)*(-85.5 +- 0)*s^2+(69.75 +- 4.749866019E-47)*s+(-137.0833333 +- 4.293552385E-46)*s^2+(6.75 +- 0)*s^(-1)+log(s)*(-18 +- 0)+log(s)*s^(-2)*(-0.5 +- 0)+log(s)*(4.5 +- 0)*s^(-1))*ep + (log(s)*s^2*(-128.25 +- 0)+(6.75 +- 0)*log(s)*s^(-1)+(7.875 +- 0)*s^(-1)+(-42.75 +- 0)*log(s)^2*s^2+(69.75 +- 0)*log(s)^2*s^3+(-9 +- 0)*log(s)^2+log(s)*(-27 +- 0)+s^3*(255.4166667 +- 1.493437348E-45)+(-154.125 +- 3.649566497E-46)*s^2+(22.5 +- 0)*log(s)^2*s+(-0.25 +- 0)*log(s)^2*s^(-2)+(-31.625 +- 0)+(209.25 +- 0)*log(s)*s^3+(-0.75 +- 0)*log(s)*s^(-2)+(-0.875 +- 0)*s^(-2)+s*(79.875 +- 4.032053447E-47)+(67.5 +- 0)*log(s)*s+(2.25 +- 0)*log(s)^2*s^(-1))*ep^2" << endl;
        
    return 0;
}
