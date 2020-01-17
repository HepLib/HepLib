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
        
//    fp.LoopMomenta = lst{q1,q2};
//    fp.Propagators = lst{
//-2*p3*q2 - pow(q2,2),
//2*p3*q2 - pow(q2,2),
//2*p1*q1 - pow(m,2) - pow(q1,2),
//2*p1*q1 + 2*p3*q1 + 2*p1*q2 + 2*p3*q2 - 2*q1*q2 - s/4 + pow(m,2) - pow(q1,2) - pow(q2,2),
//4*p1*q1 + 4*p3*q1 - s + pow(m,2) - pow(q1,2),
//4*p1*q1 + 2*p3*q1 + 4*p1*q2 + 2*p3*q2 - 2*q1*q2 - s/2 - pow(m,2) - pow(q1,2) - pow(q2,2)
//    };

    
    fp.LoopMomenta = lst{q1,q2};
    fp.Propagators = lst {
-2*p3*q2 - pow(q2,2),
pow(m,2) - pow(q1,2),
2*p3*q2 - pow(q2,2),
2*p1*q1 + 2*p3*q1 + 2*p1*q2 + 2*p3*q2 - 2*q1*q2 - s/4 + pow(m,2) - pow(q1,2) - pow(q2,2),
4*p1*q1 + 4*p3*q1 - s + pow(m,2) - pow(q1,2),
4*p1*q1 + 2*p3*q1 + 4*p1*q2 + 2*p3*q2 - 2*q1*q2 - s/2 - pow(m,2) - pow(q1,2) -pow(q2,2)
    };

    fp.Exponents = lst{ 1, 1, 1, 1, 1, 2 };
    fp.lReplacements[p1*p3] = s/8-m2;
    fp.lReplacements[p3*p3] = m2;
    fp.lReplacements[p1*p1] = m2;
    fp.lReplacements[m] = 1;
    fp.lReplacements[s] = 39;
    
    SD work;
    work.epN = -1;
    work.Verbose = 100;
    work.use_exp = true;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    work.Evaluate(fp);

    cout << work.ResultError << endl;
    
    
    return 0;
}
