#include "SD.h"
#include "mpreal.h"

using namespace HepLib;

int main(int argc, char** argv) {

    SD::debug = true;

    auto ep = SD::ep;
    
    XIntegrand xint;
    xint.Functions = lst{ 1, -x(1)-PL(0)*x(2)};
    xint.Exponents = lst{ 1, -1+ep };
    SD work;
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    work.Parameter[0] = numeric("1/3");
    work.Verbose = 100;
    work.epN = 1;
    work.CheckEnd = true;
    work.Evaluate(xint);
    cout << work.VEResult() << endl;
    
    return 0;
}
