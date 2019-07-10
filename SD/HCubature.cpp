#include "SD.h"

namespace HepLib {

/*********************************************************/
// HCubature Classes
/*********************************************************/
#include "HCubature.h"
bool HCubature::useQ(unsigned xdim, qREAL const *x) {
    if(xdim<2) return true;
    for(int i=0; i<xdim; i++) {
        if(x[i] < 1.0E-4Q) return true;
    }
    return false;
}

int HCubature::Wrapper(unsigned int xdim, long long npts, const qREAL *x, void *fdata, unsigned int ydim, qREAL *y) {
    auto self = (HCubature*)fdata;
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i=0; i<npts; i++) {
        if(useQ(xdim, x+i*xdim)) self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
        else {
            self->Integrand(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            bool ok = true;
            for(int j=0; j<ydim; j++) {
                qREAL ytmp = *(y+i*ydim+j);
                if(isnanq(ytmp) || fabsq(ytmp) > 1E8 ) {
                    ok = false;
                    break;
                }
            }
            if(!ok) self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
        }
        if(self->ReIm == 1) {
            y[i*ydim+1] = 0;
        } else if(self->ReIm == 2) {
            y[i*ydim+0] = 0;
        }
    }
    return 0;
}

void HCubature::DefaultPrintHooker(qREAL* result, qREAL* epsabs, long long int* nrun, void *fdata) {
    auto self = (HCubature*)fdata;
    
    if(self->Verbose>3 && self->RunMAX>0) {
        auto co = VE(CppFormat::q2ex(result[0]),CppFormat::q2ex(epsabs[0]))+I*VE(CppFormat::q2ex(result[1]),CppFormat::q2ex(epsabs[1]));
        cout << "     N: " << (*nrun) << ", " << VEResult(VESimplify(co)) << endl;
    }
    
    if((isnanq(result[0]) || isnanq(result[1]) || isnanq(epsabs[0]) || isnanq(epsabs[1])) && (*nrun)>self->MaxPTS/100) {
         *nrun = self->MaxPTS + 1000;
         if(self->LastState>0) self->LastState = -1;
         return;
    }
    
    if(self->RunMAX>0 && (epsabs[0] > 1E15*self->EpsAbs || epsabs[1] > 1E15*self->EpsAbs) && (*nrun)>self->MaxPTS/10) {
         *nrun = self->MaxPTS + 1000;
         if(self->LastState>0) self->LastState = -1;
         return;
    }
    
    if((self->LastState == 0) || (*nrun)<5*self->RunPTS || (epsabs[0]<=10*self->LastAbsErr[0] && epsabs[1]<=10*self->LastAbsErr[1])) {
        self->LastResult[0] = result[0];
        self->LastResult[1] = result[1];
        self->LastAbsErr[0] = epsabs[0];
        self->LastAbsErr[1] = epsabs[1];
        self->LastState = 1;
    }
    
    bool rExit = (epsabs[0] < self->EpsAbs) || (epsabs[0] < fabsq(result[0])*self->EpsRel);
    bool iExit = (epsabs[1] < self->EpsAbs+1E-50Q) || (epsabs[1] < fabsq(result[1])*self->EpsRel+1E-50Q);
    if(rExit && iExit) *nrun = self->MaxPTS + 1000;
}

ex HCubature::Integrate(unsigned int xdim, SD_Type fp, SD_Type fpQ, const qREAL* pl, const qREAL* la) {
    Integrand = fp;
    IntegrandQ = fpQ;
    Parameter = pl;
    Lambda = la;
    
    unsigned int ydim = 2;
    qREAL result[ydim], estabs[ydim];

    qREAL xmin[xdim], xmax[xdim];
    for(int i=0; i<xdim; i++) {
        xmin[i] = 0.0Q;
        xmax[i] = 1.0Q;
    }
    LastState = 0;
    
    long long run_pts = RunPTS;
    MaxPTS = RunPTS * RunMAX;
    if(MaxPTS<0) MaxPTS = -MaxPTS;
    if(xdim<2) RunPTS = 5000;
    
    hcubature_v(ydim, Wrapper, this, xdim, xmin, xmax, RunPTS, MaxPTS, EpsAbs, EpsRel, ERROR_INDIVIDUAL, result, estabs, PrintHooker);
    
    RunPTS = run_pts;
    
    if(LastState<0 && (xdim<2 || use_last)) {
        result[0] = LastResult[0];
        result[1] = LastResult[1];
        estabs[0] = LastAbsErr[0];
        estabs[1] = LastAbsErr[1];
    }
    
    ex FResult = 0;
    if(isnanq(result[0]) || isnanq(result[1])) FResult += SD::NaN;
    else {
        try{
            FResult += VE(CppFormat::q2ex(result[0]), CppFormat::q2ex(estabs[0]));
            FResult += VE(CppFormat::q2ex(result[1]), CppFormat::q2ex(estabs[1])) * I;
        } catch(...) {
            FResult += SD::NaN;
        }
    }
    
    return FResult;
}

}
