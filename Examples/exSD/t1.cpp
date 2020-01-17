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


    fp.LoopMomenta = lst {k, r, q};
    fp.Propagators = lst{ -pow(k,2),-pow(k + p1 + p2,2),-pow(-k + r,2),-pow(p1 + r,2),-pow(k - q,2),-pow(p1 + q,2) };
    fp.Exponents = lst{ 1,1,1,1,1,1 };
    fp.lReplacements[p1*p1] = 0;
    fp.lReplacements[p2*p2] = 0;
    fp.lReplacements[p2*p1] = s/2;
    fp.lReplacements[s] = 1;
    fp.Prefactor = pow(I*pow(Pi,2-ep), -3) * pow(tgamma(1-ep), 3);

    
    fp.nReplacements[ep] = ex(1)/11;
    
    SD work;
    work.epN = 0;
    work.Verbose = 2;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    work.Evaluate(fp);

    cout << work.VEResult() << endl;
    cout << "check with:" << endl;
    cout << "I*(122.72686892322528 +- 2.76707E-4)+ep^(-1)*((3.152120825132521 +- 2.93597E-5)+I*(25.132741218147014 +- 1.01913E-6))+ep^(-3)*(0.3333333333330991 +- 9.9034E-10)+(-29.193941228099536 +- 6.61155E-5)+(I*(3.141592653587586 +- 9.3337E-9)+(2.6666666655450135 +- 1.081335E-7))*ep^(-2)" << endl;
    
    return 0;
}
