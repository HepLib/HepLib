#include "SD.h"

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/FunctionMinimum.h"
#include <string>
#include <iostream>

namespace HepLib {

class FCN : public ROOT::Minuit2::FCNBase {
private:
    MinimizeBase::FunctionType InnerFunction;
    dREAL *PL;
    dREAL *LAS;
    
public:
    double operator()(const std::vector<double>& vec) const {
        int nx = vec.size();
        dREAL x[nx];
        for(int i=0; i<nx; i++) x[i] = vec[i];
        return InnerFunction(nx, x, PL, LAS);
    }

    double Up() const {
        return 1.;
    }
    
    FCN(MinimizeBase::FunctionType ff, dREAL* pl, dREAL *las) {
        InnerFunction = ff;
        PL = pl;
        LAS = las;
    }
    
};

dREAL MinUit::FindMinimum(int nvars, FunctionType func, dREAL *PL, dREAL *LAS, dREAL *UB, dREAL *LB, bool compare0) {
    double ub[nvars], lb[nvars];
    
    if(UB != NULL) for(int i=0; i<nvars; i++) ub[i] = UB[i];
    else for(int i=0; i<nvars; i++) ub[i] = 1;
    
    if(LB != NULL) for(int i=0; i<nvars; i++) lb[i] = LB[i];
    else for(int i=0; i<nvars; i++) lb[i] = 0;
    
    double step[nvars];
    for(int i=0; i<nvars; i++) step[i] = 1E-3 * (ub[i]-lb[i]);
    
    #define savePTS 5
    #define tryPTS 5
    double mPoints[savePTS][nvars], mValue[savePTS];
    for(int i=0; i<savePTS; i++) mValue[i] = 1E5;
    int max_index = 0;
    dREAL pts[tryPTS+1];
    for(int i=0; i<=tryPTS; i++) pts[i] = i*1.0/tryPTS;
    for(long long ii=0; ii<std::pow(tryPTS+1, nvars); ii++) {
        dREAL iPoints[nvars];
        int li = ii;
        for(int i=0; i<nvars; i++) {
            int mi = li % (1+tryPTS);
            iPoints[i] = lb[i] + pts[mi] * (ub[i]-lb[i]);
            li /= (1+tryPTS);
        }
        auto tmp = func(nvars, iPoints, PL, LAS);
        if(compare0 && tmp<-1E-50) return -1;
        if(mValue[max_index] > tmp) {
            mValue[max_index] = tmp;
            for(int j=0; j<nvars; j++) mPoints[max_index][j] = iPoints[j];
            max_index = 0;
            for(int j=0; j<savePTS; j++) {
                if(mValue[j] > mValue[max_index]) max_index = j;
            }
        }
    }
    
    double ret = 1E5;
    for(int ii=0; ii<savePTS; ii++) {
        FCN fcn(func, PL, LAS);
        ROOT::Minuit2::MnUserParameters upar;
        for(int i=0; i<nvars; i++) {
            ostringstream xs;
            xs << "x" << i;
            upar.Add(xs.str(), mPoints[ii][i], step[i], lb[i], ub[i]);
        }
 
        ROOT::Minuit2::MnMigrad minizer(fcn, upar);
        ROOT::Minuit2::FunctionMinimum fmin = minizer();
        auto tmp_ret = fmin.Fval();
        if(tmp_ret < ret) ret = tmp_ret;
        if(compare0 && ret<-1E-50) return -1;
    }
    return ret;
}

}
