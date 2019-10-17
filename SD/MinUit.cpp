#include "SD.h"

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/FunctionMinimum.h"
#include <string>
#include <iostream>
#include <cstdlib>
#include <ctime>

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

int MinUit::TryPTS = 5;
int MinUit::SavePTS = 2;

void MinUit_Random(int n, double *x) {
    static bool inited = false;
    if(!inited) {
        srand(static_cast<unsigned int>(clock()));
    }
    for(int i=0; i<n; i++) {
        x[i] = double(rand())/(double(RAND_MAX)+1.0);
    }
}

dREAL MinUit::FindMinimum(int nvars, FunctionType func, dREAL *PL, dREAL *LAS, dREAL *UB, dREAL *LB, dREAL *IP, bool compare0) {
    double ub[nvars], lb[nvars];
    
    if(UB != NULL) for(int i=0; i<nvars; i++) ub[i] = UB[i];
    else for(int i=0; i<nvars; i++) ub[i] = 1;
    
    if(LB != NULL) for(int i=0; i<nvars; i++) lb[i] = LB[i];
    else for(int i=0; i<nvars; i++) lb[i] = 0;
    
    double step[nvars];
    for(int i=0; i<nvars; i++) step[i] = 5E-4 * (ub[i]-lb[i]);
    
    double minPoints[SavePTS][nvars], maxPoints[SavePTS][nvars], minValue[SavePTS], maxValue[SavePTS];
    for(int i=0; i<SavePTS; i++) {
        minValue[i] = 1E5;
        maxValue[i] = -1E5;
    }
    int min_index=0, max_index=0;
    dREAL pts[TryPTS+1];
    for(int i=0; i<=TryPTS; i++) pts[i] = i*1.0/TryPTS;
    for(long long ii=0; ii<std::pow(TryPTS+1, nvars); ii++) {
        dREAL iPoints[nvars];
        int li = ii;
        for(int i=0; i<nvars; i++) {
            int mi = li % (1+TryPTS);
            iPoints[i] = lb[i] + pts[mi] * (ub[i]-lb[i]);
            li /= (1+TryPTS);
        }
        auto tmp = func(nvars, iPoints, PL, LAS);
        if(compare0 && tmp<ZeroValue) return -1;
        
        if(minValue[max_index] > tmp) {
            minValue[max_index] = tmp;
            for(int j=0; j<nvars; j++) minPoints[max_index][j] = iPoints[j];
            max_index = 0;
            for(int j=0; j<SavePTS && j<std::pow(TryPTS+1, nvars); j++) {
                if(minValue[j] > minValue[max_index]) max_index = j;
            }
        }
        
        if(maxValue[min_index] < tmp) {
            maxValue[min_index] = tmp;
            for(int j=0; j<nvars; j++) maxPoints[max_index][j] = iPoints[j];
            min_index = 0;
            for(int j=0; j<SavePTS && j<std::pow(TryPTS+1, nvars); j++) {
                if(maxValue[j] < maxValue[max_index]) min_index = j;
            }
        }
        
    }
    
    double ret = 1E5;
    
    for(int ii=0; ii<2*SavePTS; ii++) {
        FCN fcn(func, PL, LAS);
        ROOT::Minuit2::MnUserParameters upar;
        double rndPoints[nvars];
        MinUit_Random(nvars, rndPoints);
        for(int i=0; i<nvars; i++) {
            ostringstream xs;
            xs << "x" << i;
            upar.Add(xs.str(), rndPoints[i], step[i], lb[i], ub[i]);
        }
 
        ROOT::Minuit2::MnMigrad minizer(fcn, upar, 2);
        ROOT::Minuit2::FunctionMinimum fmin = minizer();
        auto tmp_ret = fmin.Fval();
        if(tmp_ret < ret) ret = tmp_ret;
        if(compare0 && ret<ZeroValue) return -1;
    }
    
    for(int ii=0; ii<SavePTS && ii<std::pow(TryPTS+1, nvars); ii++) {
        FCN fcn(func, PL, LAS);
        ROOT::Minuit2::MnUserParameters upar;
        for(int i=0; i<nvars; i++) {
            ostringstream xs;
            xs << "x" << i;
            upar.Add(xs.str(), minPoints[ii][i], step[i], lb[i], ub[i]);
        }
 
        ROOT::Minuit2::MnMigrad minizer(fcn, upar, 2);
        ROOT::Minuit2::FunctionMinimum fmin = minizer();
        auto tmp_ret = fmin.Fval();
        if(tmp_ret < ret) ret = tmp_ret;
        if(compare0 && ret<ZeroValue) return -1;
    }
    
    for(int ii=0; ii<SavePTS && ii<std::pow(TryPTS+1, nvars); ii++) {
        FCN fcn(func, PL, LAS);
        ROOT::Minuit2::MnUserParameters upar;
        for(int i=0; i<nvars; i++) {
            ostringstream xs;
            xs << "x" << i;
            upar.Add(xs.str(), maxPoints[ii][i], step[i], lb[i], ub[i]);
        }
 
        ROOT::Minuit2::MnMigrad minizer(fcn, upar);
        ROOT::Minuit2::FunctionMinimum fmin = minizer();
        auto tmp_ret = fmin.Fval();
        if(tmp_ret < ret) ret = tmp_ret;
        if(compare0 && ret<ZeroValue) return -1;
    }
    
    return ret;
}

void MinUit::ForceStop() {
    throw domain_error("force stop!");
}

void MinUit::Minimize(int nvars, FunctionType func, dREAL *ip) {
    try {
        omp_set_num_threads(1);
        FCN fcn(func, NULL, NULL);
        ROOT::Minuit2::MnUserParameters upar;
        for(int i=0; i<nvars; i++) {
            ostringstream xs;
            xs << "x" << i;
            upar.Add(xs.str(), ip[i], 0.5, 0, 10);
        }

        ROOT::Minuit2::MnMigrad minizer(fcn, upar);
        minizer();
    } catch(...) { }
}

}
