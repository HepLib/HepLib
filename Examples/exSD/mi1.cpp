#include "SD.h"

using namespace HepLib;

int main(int argc, char** argv) {

    //SD::debug = true;
    
    auto ep = SD::ep;
    auto iep = SD::iEpsilon;
    
    Digits = 30;
    
    symbol k1("k1"), k2("k2"), q1("q1"), q2("q2");
    ex mc = 1;
    
    FeynmanParameter fp;
    
    fp.LoopMomenta = lst {q1, q2};
    fp.Propagators = lst{-pow(q1, 2), 2*k2*q1 - pow(q1, 2), k1*q1 + k2*q1 - pow(q1, 2),
        -(k1*q1) + k2*q1 + k1*q2 - k2*q2 + 2*q1*q2 + pow(mc, 2) - pow(q1, 2)
        - pow(q2, 2), k1*q1 + k2*q1 - k1*q2 - k2*q2 + 2*q1*q2 - pow(mc, 2) -
        pow(q1, 2) - pow(q2, 2), pow(mc, 2) - pow(q2, 2), k1*q2 + k2*q2 -
        pow(mc, 2) - pow(q2, 2)};
    fp.Exponents = lst{ 1, 1, 1, 1, 1, 1, 1 };
    fp.lReplacements[k1*k1] = 0;
    fp.lReplacements[k2*k2] = 0;
    fp.lReplacements[k1*k2] = 2*mc*mc;
    //fp.Prefactor = pow(I*pow(Pi,2-ep)*exp(-ep*Euler), -1);
        
    SD work;
    work.epN = 0;
    work.Verbose = 100;
    //work.ParallelProcess = 0;
    
    //work.use_cpp = false;
    //SD::debug = true;
    
    char *CFLAGS = getenv("SD_CFLAGS");
    work.CFLAGS = CFLAGS;
    
    cout << endl << "Starting @ " << now() << endl;
    if(work.SecDec==NULL) work.SecDec = new SecDecG();
    if(work.Integrator==NULL) work.Integrator = new HCubature();
    if(work.Minimizer==NULL) work.Minimizer = new MinUit();
    if(strlen(CFLAGS)<1) CFLAGS = getenv("SD_CFLAGS");
    
    work.Initialize(fp);
    if(work.FunExp.size()<1) return 0;
    work.Scalelesses();
    work.ChengWu();
for(auto kv : work.FunExp) {
    cout << kv.first.op(1) << endl;
}
exit(0);
    work.RemoveDeltas();
    work.KillPowers();
    work.SDPrepares();
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
