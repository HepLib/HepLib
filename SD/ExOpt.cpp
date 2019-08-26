#include "SD.h"

#include "optim.hpp"

namespace HepLib {

double ExOpt_ObjectWrapper(const arma::vec& x, arma::vec* grad_out, void* vself) {
    ExOpt *self = (ExOpt*)vself;
    int nvars = x.size();
    dREAL xa[nvars];
    for(int i=0; i<nvars; i++) xa[i] = x(i);
    double ret = self->ObjectFunction(nvars, xa, self->PL, self->LAS);
    if(self->Compare0 && ret < self->ZeroValue) self->ForceStop();
    return ret;
}

void ExOpt::ForceStop() {
    throw domain_error("force stop!");
}

dREAL ExOpt::FindMinimum(int nvars, FunctionType func, dREAL *pl, dREAL *las, dREAL *ub, dREAL *lb, dREAL *IP, bool compare0) {
    ObjectFunction = func;
    PL = pl;
    LAS = las;
    Compare0 = compare0;
    
    arma::vec x0(nvars), slb(nvars), sub(nvars);
    if(ub != NULL) for(int i=0; i<nvars; i++) sub(i) = ub[i];
    else for(int i=0; i<nvars; i++) sub(i) = 1;
    if(lb != NULL) for(int i=0; i<nvars; i++) slb(i) = lb[i];
    else for(int i=0; i<nvars; i++) slb(i) = 0;
    for(int i=0; i<nvars; i++) x0(i) = (sub(i)+slb(i))/2.0;
    
    optim::algo_settings_t settings;
    settings.vals_bound = true;
    settings.lower_bounds = slb;
    settings.upper_bounds = sub;
    settings.err_tol = 1E-3;
    settings.iter_max = 100;
    settings.de_n_gen = 30;
    settings.de_max_fn_eval = 500;
    
    dREAL ret;
    
    try {
        // de pso de_prmm
        bool success = optim::de(x0, ExOpt_ObjectWrapper, this, settings);
        if(success) {
            ret = ExOpt_ObjectWrapper(x0, NULL, this);
        } else {
            cout << "Failed ..." << endl;
        }
    } catch(domain_error) {
        return -1;
    }
    cout << ret << endl;
    return ret;
}

void ExOpt::Minimize(int nvars, FunctionType func, dREAL *ip) {
    ObjectFunction = func;
    
    arma::vec x0(nvars), slb(nvars), sub(nvars);
    for(int i=0; i<nvars; i++) {
        sub(i) = 10;
        slb(i) = 1E-4;
        x0(i) = ip[i];
    };
    
    optim::algo_settings_t settings;
    settings.vals_bound = true;
    settings.lower_bounds = slb;
    settings.upper_bounds = sub;
    
    settings.err_tol = 1E-3;
    settings.iter_max = 100;
    settings.de_n_gen = nvars;
    settings.de_max_fn_eval = 100;
    
    dREAL ret;
    
    try {
        optim::de(x0, ExOpt_ObjectWrapper, this, settings);
    } catch(domain_error) { }
}


}
