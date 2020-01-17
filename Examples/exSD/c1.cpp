#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    symbol k("k"), l("l"), r("r"), q("q"), q1("q1"), q2("q2");
    symbol p1("p1"), p2("p2"), p3("p3"), p4("p4"), p5("p5");
    symbol m("m"), s("s"), t("t"), u("u");
    
    ex m2=m*m;
    XIntegrand xint;
    
    xint.Functions = lst{ 1, x(0)*x(2)-3 };
    xint.Exponents = lst{ 1, -1+ep };
    
    SD work;
    work.epN = 1;
    work.Verbose = 2;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    work.Evaluate(xint);

    cout << work.VEResult() << endl;

    return 0;
}
