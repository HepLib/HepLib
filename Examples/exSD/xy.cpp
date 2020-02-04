#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
            
    SD work;
    work.epN = 0;
    work.Verbose = 100;
    work.CheckEnd = true;
    
    cout << endl << "Starting @ " << now() << endl;
    if(work.SecDec==NULL) work.SecDec = new SecDecG();
    if(work.Integrator==NULL) work.Integrator = new HCubature();
    if(work.Minimizer==NULL) work.Minimizer = new MinUit();
    
    work.FunExp.clear();
    lst fun = lst{ -x(1)/2+ex(1)/2, -x(1)/2 };
    lst exp = lst{ 1,  -2+ep };
    
    work.FunExp.push_back(make_pair(fun, exp));
    //work.Deltas.push_back(lst{x(0),x(1),x(2)});

//work.Normalizes();
//work.RemoveDeltas();
////work.KillPowers();
//work.Normalizes();
//work.SDPrepares();
//for(auto kv : work.Integrands) {
//    cout << kv << endl;
//}
//exit(0);
    
    //work.FunExp = SD::Binarize(make_pair(fun, exp), x(2)-16*x(1));
        
    
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    
    
    
    work.KillPowers();
    work.ChengWu();
    //exit(0);

    work.Normalizes();
    work.RemoveDeltas();
    //work.KillPowers();
    work.Normalizes();

    work.SDPrepares();
//exit(0);
    work.EpsEpExpands();
    work.CIPrepares();
    auto pps = work.ParallelProcess;
    work.ParallelProcess = 0;
    work.Contours();
    work.Integrates();
    work.ParallelProcess = pps;
    delete work.SecDec;
    delete work.Integrator;
    delete work.Minimizer;
    cout << "Finished @ " << now() << endl << endl;

    cout << work.VEResult() << endl;
        
    return 0;
}
