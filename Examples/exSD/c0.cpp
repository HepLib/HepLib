#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {

    //SD::debug = true;
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    Digits = 30;
    
    symbol k("k"), l("l"), r("r"), q("q"), q1("q1"), q2("q2");
    symbol p1("p1"), p2("p2"), p3("p3"), p4("p4"), p5("p5");
    symbol m("m"), s("s"), t("t"), u("u");
    
    ex m2=m*m;
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
    fp.lReplacements[s] = -3;
    fp.lReplacements[t] = 1;
    fp.Prefactor = pow(I*pow(Pi,2-ep)*exp(-ep*Euler), -1);
        
    SD work;
    work.epN = 0;
    work.Verbose = 100;
    //work.ParallelProcess = 0;
    
    //work.use_cpp = false;
    //SD::debug = true;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    if(false) {
        auto cuba = new CUBA();
        cuba->Method = CUBA::CUHRE;
        cuba->VEGAS_NSTART = 1000;
        cuba->VEGAS_NINCREASE = 1000;
        cuba->VEGAS_NBATCH = 1000;
        work.Integrator = cuba;
        work.TryPTS = 1000000;
        work.RunPTS = 1000000;
        work.RunMAX = 10;
    }
    work.Evaluate(fp);

    cout << work.VEResult() << endl;
    cout << "check with:" << endl;
    cout << "I*(-6.6265E-11 +- 1.61138E-6)+(4.386490845011835 +- 1.61436E-6)+(I*(-8.821394565321117E-20 +- 6.43574E-24)+(-1.33333333333333333339295812540587255201 +- 5.5055E-24))*ep^(-2)+ep^(-1)*(I*(-2.09439510241206 +- 3.77364E-7)+(0.7324081924598961 +- 3.73657E-7))" << endl;
        
    return 0;
}
