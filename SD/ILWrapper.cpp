#include "SD.h"

namespace HepLib {

int ILWrapper::Verbose;
unsigned int ILWrapper::xsize;
IntegratorBase::SD_Type ILWrapper::fp = NULL;
IntegratorBase::SD_Type ILWrapper::fpQ = NULL;
IntegratorBase *ILWrapper::Integrator = NULL;
MinimizeBase *ILWrapper::miner = NULL;
qREAL *ILWrapper::paras = NULL;
dREAL *ILWrapper::lambda = NULL;
dREAL ILWrapper::err_max;
dREAL ILWrapper::hjRHO = 0.75;
dREAL ILWrapper::err_min = -0.001;
long long ILWrapper::MaxPTS = 500;
long long ILWrapper::RunPTS = 0;
ex ILWrapper::lastResErr;

dREAL ILWrapper::IntError(int nvars, dREAL *las, dREAL *n1, dREAL *n2) {
    RunPTS++;
    qREAL qlas[nvars];
    for(int i=0; i<nvars; i++) qlas[i] = las[i];
    auto res = Integrator->Integrate(xsize, fp, fpQ, paras, qlas);
    if(res.has(SD::NaN)) return 1.E15;
    auto err = res.subs(VE(wild(0), wild(1))==wild(1));
    numeric nerr = numeric("1E15");
    try {
        nerr = ex_to<numeric>(abs(err).evalf());
        if(nerr > numeric("1E15")) nerr = numeric("1E15");
        if(err_max > nerr.to_double()) {
            auto diff = VESimplify(lastResErr - res);
            diff = diff.subs(VE(0,0)==0);
            exset ves;
            diff.find(VE(wild(0), wild(1)), ves);
            for(auto ve : ves) {
                if(abs(ve.op(0)) > ve.op(1)) return 1E15;
            }
            if(Verbose > 3) {
                cout << "\r                             \r";
                cout << WHITE << "     " << RunPTS << ": " << RESET;
                for(int i=0; i<nvars; i++) cout << las[i] << " ";
                cout << endl << "     " << res.subs(VE(0,0)==0) << endl;
            }
            err_max = nerr.to_double();
            for(int i=0; i<nvars; i++) lambda[i] = las[i];
            if(err_max<=err_min) {
                cout << "\r                             \r";
                cout << "     ------------------------------" << endl;
                miner->ForceStop();
            }
        } else {
            if(Verbose > 3) {
                cout << "\r                             \r";
                cout << WHITE << "     [ " << RunPTS << " / " << MaxPTS << " ] ..." << RESET << flush;
            }
        }
    } catch(domain_error &ex) {
        throw ex;
    } catch(...) { }
    if(RunPTS>=MaxPTS) {
        cout << "\r                             \r";
        cout << "     ------------------------------" << endl;
        miner->ForceStop();
    }
    return nerr.to_double();
}


}
