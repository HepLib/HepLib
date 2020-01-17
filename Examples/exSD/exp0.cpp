#include "SD.h"
#include <cmath>

using namespace HepLib;

int main(int argc, char** argv) {

    SD::debug = true;
    
    auto ep = SD::ep;
    auto eps = SD::eps;
    auto iep = SD::iEpsilon;
    auto t = SD::vs;

    symbol k("k"), l("l"), r("r"), q("q"), q1("q1"), q2("q2");
    symbol p1("p1"), p2("p2"), p3("p3"), p4("p4"), p5("p5");
    symbol m("m"), s("s"), u("u");

    
    Digits = 30;
        
    FeynmanParameter fp;
    
    fp.LoopMomenta = lst {k};
    fp.Propagators = lst{ -pow(k,2),-pow(k + p1,2),-pow(k + p1 + p2,2),-pow(k + p1 + p2 + p4,2) };
    fp.Exponents = lst{ 1, 1, 1, 1 };
    fp.lReplacements[p1*p1] = 0;
    fp.lReplacements[p2*p2] = 0;
    fp.lReplacements[p4*p4] = 0;
    fp.lReplacements[p1*p2] = -s/2;
    fp.lReplacements[p2*p4] = -t/2;
    fp.lReplacements[p1*p4] = (s+t)/2;
    fp.lReplacements[s] = 1;
    fp.Prefactor = pow(I*pow(Pi,2-ep)*exp(-ep*Euler), -1);
    
    fp.nReplacements[ep] = ex(1)/11;
    
    SD work;
    work.epN = 1;
    work.sN = 1;
    work.Verbose = 100;
    //work.ParallelProcess = 0;
    work.RunPTS = 1000;
    work.RunMAX = 1000;
    
    auto hc = new HCubature();
    hc->use_last = true;
    work.Integrator = hc;
    
    //work.use_cpp = false;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    work.Evaluate(fp);

    work.VEPrint();
    
    cout << WHITE << "Check with: " << RESET << endl;
    
    cout << "(s^(-1.0)*(4.0 +- 0))*ep^(-2) + ((0 +- 1.56929E-15)+s*(0 +- 1.39681E-12)+s^(-1)*(-2.0 +- 0)*log(s))*ep^(-1) + ((0 +- 4.0313E-6)+s*(0 +- 5.4998E-6)+s^(-1.0)*(-13.15947253478581149177932135310982312 +- 4.02865E-18)) + ((1.9999999882487511 +- 6.28736E-6)*log(s)+s*(5.184801111112736 +- 1.33817E-5)+s^(-1.0)*(-13.623311572379617 +- 7.4906E-6)+(-1.0 +- 5.5483E-16)*log(s)^2+(-11.869605299256339 +- 1.76101E-5)+s^(-1)*(0.333333333333333333333333333333333333334 +- 0)*log(s)^3+s*(0.5 +- 4.9385E-13)*log(s)^2+s^(-1)*(11.51453846793758505530690618147839252925 +- 3.02147E-18)*log(s)+s*(-0.4999999936021155 +- 4.6124E-6)*log(s))*ep" << endl;
        
    return 0;
}
