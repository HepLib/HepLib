#include "SD.h"
#include "mpreal.h"

using namespace HepLib;

int main(int argc, char** argv) {

    //SD::debug = true;

    auto ep = SD::ep;
    
    if(false) {
        XIntegrand xint;
        xint.Functions = lst{1, (-(x(1)*x(2)*(-9 + 92*pow(x(3), 2) - 108*x(3))) + 8*pow(x(2), 2)*(3 + 2*pow(x(3), 2) - 5*x(3))*x(3) + 6*pow(x(1), 2)*(3 + 22*x(3)))/27, x(1), x(2), 1 - x(3), x(3), 2*x(1) + x(2) - (2*x(2)*x(3))/3, (-4*pow(x(2), 2)*(-1 + x(3))*x(3) + x(1)*x(2)*(3 + 24*pow(x(3), 2) - 32*pow(x(3), 3) + 30*x(3)) + pow(x(1), 2)*(6 + 92*pow(x(3), 2) + 52*x(3)))/27, 1};
        xint.Exponents = lst{-1, -2*(1 + ep), -2 + ep, 1, -ep, ep, -1 + 3*ep, 2, -2 - 3*ep};
        xint.Deltas.push_back(lst{x(1),x(2)});
        SD work;
        char *CFLAGS = getenv("SD_CFLAGS");
        work.CFLAGS = CFLAGS;
        //work.ParallelProcess = 0;
        work.Verbose = 0;
        work.epN = 0;
        work.CheckEnd = true;
        
        work.TryPTS = 1000000;
        work.RunPTS = 1000000;
        work.RunMAX = 10;
        
        work.Evaluate(xint);
        cout << work.VEResult() << endl;
    }
    
    
    
    return 0;
}
