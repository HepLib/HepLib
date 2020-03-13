#include "SD.h"

namespace HepLib::SD {

int ErrMin::Verbose;
IntegratorBase *ErrMin::Integrator = NULL;
MinimizeBase *ErrMin::miner = NULL;
qREAL *ErrMin::paras = NULL;
dREAL *ErrMin::lambda = NULL;
dREAL ErrMin::err_max;
dREAL ErrMin::hjRHO = 0.75;
dREAL ErrMin::err_min = -0.001;
long long ErrMin::MaxRND = 50;
long long ErrMin::RunRND = 0;
ex ErrMin::lastResErr;

dREAL ErrMin::IntError(int nvars, dREAL *las, dREAL *n1, dREAL *n2) {
    RunRND++;
    auto digits = Digits;
    Digits = 50;
    qREAL qlas[nvars];
    for(int i=0; i<nvars; i++) qlas[i] = las[i];
    Integrator->Lambda = qlas;
    auto res = Integrator->Integrate();
    if(res.has(NaN)) {
        Digits = digits;
        return 1.E100;
    }
    auto err = res.subs(VE(w0, w1)==w1);
    numeric nerr = numeric(1.E100);
    try {
        nerr = ex_to<numeric>(abs(err).evalf());
        if(nerr > numeric(1.E100)) nerr = numeric(1.E100);
        if(err_max > nerr.to_double()) {
            auto diff = VESimplify(lastResErr - res);
            diff = diff.subs(VE(0,0)==0);
            exset ves;
            diff.find(VE(w0, w1), ves);
            for(auto ve : ves) {
                if(abs(ve.op(0)) > ve.op(1)) {
                    Digits = digits;
                    return 1.E100;
                }
            }
            if(Verbose > 3) {
                cout << "\r                             \r";
                cout << Color_HighLight << "     " << RunRND << ": " << RESET;
                for(int i=0; i<nvars; i++) cout << las[i] << " ";
                cout << endl << "     " << res.subs(VE(0,0)==0).subs(VE(w1,w2)==VEO(w1,w2)) << endl;
            }
            err_max = nerr.to_double();
            for(int i=0; i<nvars; i++) lambda[i] = las[i];
            if(err_max<=err_min) {
                cout << "\r                             \r";
                cout << "     ------------------------------" << endl;
                Digits = digits;
                miner->ForceStop();
                return 0.;
            }
        } else {
            if(Verbose > 3) {
                cout << "\r                             \r";
                cout << Color_HighLight << "     [ " << RunRND << " / " << MaxRND << " ] ..." << RESET << flush;
            }
        }
    } catch(domain_error &ex) {
        throw ex;
    } catch(...) { }
    if(RunRND>=MaxRND) {
        cout << "\r                             \r";
        cout << "     ------------------------------" << endl;
        Digits = digits;
        miner->ForceStop();
        return 0.;
    }
    
    auto pid = getpid();
    ostringstream fn;
    fn << pid << ".las.done";
    if(file_exists(fn.str().c_str())) {
        ostringstream cmd;
        cmd << "rm " << fn.str();
        system(cmd.str().c_str());
        cout << "\r                             \r";
        if(Verbose>3) cout << "     Exit: " << fn.str() << endl;
        cout << "     ------------------------------" << endl;
        Digits = digits;
        miner->ForceStop();
        return 0.;
    }
    dREAL ret = nerr.to_double();
    Digits = digits;
    return ret;
}


}
