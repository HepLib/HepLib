#include "HepLib.h"

using namespace HepLib;
using namespace SD;

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
    fp.Propagator = lst{ -pow(k,2),-pow(k + p1 + p2,2),-pow(-k + r,2),-pow(p1 + r,2),-pow(k - q,2),-pow(p1 + q,2) };
    fp.Exponent = lst{ 1+3*ep,1,1,1,1,1 };
    fp.lReplacement[p1*p1] = 0;
    fp.lReplacement[p2*p2] = 0;
    fp.lReplacement[p2*p1] = s/2;
    fp.lReplacement[s] = -1;
    fp.Prefactor = pow(I*pow(Pi,2-ep), -3) * pow(tgamma(1-ep), 3);

    SecDec work;
    work.eps_lst = lst{ lst{ep, 0} };
    Verbose = 2;

    work.Evaluate(fp);

    cout << work.VEResult() << endl;
    cout << "check with:" << endl;
    cout << "1.812319(2)E1*ep^(-1)+1.2532095(6)E2+1.666666667(4)E-1*ep^(-3)+1.83333333(6)*ep^(-2)" << endl;
    
    
    return 0;
}
