#include "SD.h"

namespace HepLib {

/*********************************************************/
// HCubature Classes
/*********************************************************/
#include "HCubature.h"
bool HCubature::useQ(unsigned xdim, qREAL const *x) {
    if(xdim<2) return true;
    for(int i=0; i<xdim; i++) {
        if(x[i] < 5.0E-6Q) return true;
    }
    return false;
}

int HCubature::Wrapper(unsigned int xdim, long long npts, const qREAL *x, void *fdata, unsigned int ydim, qREAL *y) {
    auto self = (HCubature*)fdata;
    bool NaNQ = false;
    #pragma omp parallel for num_threads(omp_get_num_procs()) schedule(dynamic, 1)
    for(int i=0; i<npts; i++) {
        qREAL xx[xdim], det = 1.Q;
        for(int xi=0; xi<xdim; xi++) {
            xx[xi] = powq(x[xi+i*xdim], self->XN);
            det *= self->XN * powq(x[xi+i*xdim], self->XN-1);
        }
        if(self->UseQ || useQ(xdim, xx)) {
            self->IntegrandQ(xdim, xx, ydim, y+i*ydim, self->Parameter, self->Lambda);
            //self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
        } else {
            self->Integrand(xdim, xx, ydim, y+i*ydim, self->Parameter, self->Lambda);
            //self->Integrand(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            bool ok = true;
            for(int j=0; j<ydim; j++) {
                qREAL ytmp = y[i*ydim+j];
                if(isnanq(ytmp)) {
                    ok = false;
                    break;
                }
            }
            if(!ok) self->IntegrandQ(xdim, xx, ydim, y+i*ydim, self->Parameter, self->Lambda);
            //if(!ok) self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
        }
        
        for(int j=0; j<ydim; j++) {
            qREAL ytmp = y[i*ydim+j];
            if(isnanq(ytmp)) {
                NaNQ = true;
                break;
            }
        }
        
        if(self->ReIm == 1) {
            y[i*ydim+1] = 0;
        } else if(self->ReIm == 2) {
            y[i*ydim+0] = 0;
        }
        
        y[i*ydim+0] *= det;
        y[i*ydim+1] *= det;
    }
    return NaNQ ? 1 : 0;
}

void HCubature::DefaultPrintHooker(qREAL* result, qREAL* epsabs, long long int* nrun, void *fdata) {
    auto self = (HCubature*)fdata;
    
    if((isnanq(result[0]) || isnanq(result[1]) || isnanq(epsabs[0]) || isnanq(epsabs[1]))) {
         *nrun = self->MaxPTS + 1000;
         if(self->LastState>0) self->LastState = -1;
         return;
    }
    
    if(self->RunMAX>0 && (epsabs[0] > 1E25*self->EpsAbs || epsabs[1] > 1E25*self->EpsAbs)) {
         *nrun = self->MaxPTS + 1000;
         if(self->LastState>0) self->LastState = -1;
         return;
    }
    
    if((self->LastState == 0) || (epsabs[0]<=2*self->LastAbsErr[0] && epsabs[1]<=2*self->LastAbsErr[1])) {
        self->LastResult[0] = result[0];
        self->LastResult[1] = result[1];
        self->LastAbsErr[0] = epsabs[0];
        self->LastAbsErr[1] = epsabs[1];
        self->LastState = 1;
    }
    
    if(self->Verbose>10 && self->RunMAX>0) {
        char r0[64], r1[64], e0[32], e1[32];
        quadmath_snprintf(r0, sizeof r0, "%.10Qg", result[0]);
        quadmath_snprintf(r1, sizeof r1, "%.10Qg", result[1]);
        quadmath_snprintf(e0, sizeof e0, "%.5Qg", epsabs[0]);
        quadmath_snprintf(e1, sizeof e1, "%.5Qg", epsabs[1]);
        cout << "     N: " << (*nrun) << ", ";
        if(self->ReIm==3 || self->ReIm==1) cout << "["<<r0 << ", " << e0 << "]";
        if(self->ReIm==3 || self->ReIm==2) cout << "+I*[" << r1 << ", " << e1 << "]";
        cout << endl;
    }
    
    bool rExit = (epsabs[0] < self->EpsAbs) || (epsabs[0] < fabsq(result[0])*self->EpsRel);
    bool iExit = (epsabs[1] < self->EpsAbs+1E-50Q) || (epsabs[1] < fabsq(result[1])*self->EpsRel+1E-50Q);
    if(rExit && iExit) {
        *nrun = self->MaxPTS + 1000;
        return;
    }
    
    auto pid = getpid();
    ostringstream fn;
    fn << pid << ".int.done";
    if(file_exists(fn.str().c_str())) {
        *nrun = self->MaxPTS + 1000;
        ostringstream cmd;
        cmd << "rm " << fn.str();
        system(cmd.str().c_str());
    }
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
    
    MaxPTS = RunPTS * RunMAX;
    if(MaxPTS<0) MaxPTS = -MaxPTS;
    
    int nok = hcubature_v(ydim, Wrapper, this, xdim, xmin, xmax, RunPTS, MaxPTS, EpsAbs, EpsRel, result, estabs, PrintHooker);
    
    if(nok) {
        if( (cabsq(result[0]+result[1]*1.Qi) < FLT128_EPSILON) && (cabsq(estabs[0]+estabs[1]*1.Qi) < FLT128_EPSILON) ) {
            cout << RED << "HCubature Failed with 0 result returned!" << RESET << endl;
            return SD::NaN;
        }
    }
        
    if(LastState==-1 && use_last) {
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
