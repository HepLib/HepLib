#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    symbol k("k"), l("l"), r("r"), q("q"), q1("q1"), q2("q2");
    symbol p1("p1"), p2("p2"), p3("p3"), p4("p4"), p5("p5");
    symbol m("m"), s("s"), t("t"), u("u");
    
    ex m2=m*m;
    FeynmanParameter fp;


    //SecDec Example - 7_epsprops_triangle_3L
    fp.LoopMomenta = lst {k, r, q};
    fp.Propagators = lst{ -pow(k,2),-pow(k + p1 + p2,2),-pow(-k + r,2),-pow(p1 + r,2),-pow(k - q,2),-pow(p1 + q,2) };
    fp.Exponents = lst{ 1+3*ep,1,1,1,1,1 };
    fp.lReplacements[p1*p1] = 0;
    fp.lReplacements[p2*p2] = 0;
    fp.lReplacements[p2*p1] = s/2;
    fp.lReplacements[s] = -1;
    fp.Prefactor = pow(I*pow(Pi,2-ep), -3) * pow(tgamma(1-ep), 3);

    
    fp.nReplacements[ep] = ex(1)/11;
    
    SD work;
    work.epN = 0;
    work.Verbose = 2;
    //work.PoleRequested = -1;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    //work.Integrator = new CUBA();
    work.Evaluate(fp);

    cout << work.VEResult() << endl;
    cout << "check with:" << endl;
    cout << "(1.8333333326646433 +- 5.1808E-8)*ep^(-2)+(18.12319251247844 +- 2.0688E-5)*ep^(-1)+(0.16666666666654956 +- 4.9517E-10)*ep^(-3)+(125.3209681704151 +- 6.8991E-5)" << endl;
    
    
    return 0;
}
