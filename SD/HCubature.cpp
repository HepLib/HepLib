#include "SD.h"
#include "mpreal.h"

namespace HepLib {

/*-----------------------------------------------------*/
// HCubature Classes
/*-----------------------------------------------------*/
#include "HCubature.h"

int HCubature::Wrapper(unsigned int xdim, long long npts, const qREAL *x, void *fdata, unsigned int ydim, qREAL *y) {
    auto self = (HCubature*)fdata;
    bool NaNQ = false;

    if(self->UseCpp) {
        #pragma omp parallel for num_threads(omp_get_num_procs()) schedule(dynamic, 1)
        for(int i=0; i<npts; i++) {
            mpfr::mpreal::set_default_prec(mpfr::digits2bits(self->MPDigits));
            int iDQMP = self->inDQMP(x+i*xdim);
            if( (self->IntegrandMP!=NULL) && (self->DQMP>2 || iDQMP>2) ) {
                self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            } else if(self->DQMP>1 || iDQMP>1) {
                self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
                bool ok = true;
                for(int j=0; j<ydim; j++) {
                    qREAL ytmp = y[i*ydim+j];
                    if(isnanq(ytmp) || isinfq(ytmp)) {
                        ok = false;
                        break;
                    }
                }
                if(!ok && (self->IntegrandMP!=NULL)) self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            } else {
                self->Integrand(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
                bool ok = true;
                for(int j=0; j<ydim; j++) {
                    qREAL ytmp = y[i*ydim+j];
                    if(isnanq(ytmp) || isinfq(ytmp)) {
                        ok = false;
                        break;
                    }
                }
                if(!ok) self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
                
                ok = true;
                for(int j=0; j<ydim; j++) {
                    qREAL ytmp = y[i*ydim+j];
                    if(isnanq(ytmp) || isinfq(ytmp)) {
                        ok = false;
                        break;
                    }
                }
                if(!ok && (self->IntegrandMP!=NULL)) self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            }
            
            // Final Check NaN/Inf
            bool ok = true;
            for(int j=0; j<ydim; j++) {
                qREAL ytmp = y[i*ydim+j];
                if(isnanq(ytmp) || isinfq(ytmp)) {
                    ok = false;
                    break;
                }
            }
            if(!ok && (self->IntegrandMP!=NULL)) {
                qREAL xx[xdim];
                for(int ii=0; ii<xdim; ii++) xx[ii] = x[i*xdim+ii] < 1.Q-30 ? 1.Q-30  : x[i*xdim+ii] * 0.95Q;
                self->IntegrandMP(xdim, xx, ydim, y+i*ydim, self->Parameter, self->Lambda);
            }
        }
    } else {
        self->Integrand(xdim+npts*100, x, ydim+npts*100, y, self->Parameter, self->Lambda);
    } 
    
    for(int i=0; i<npts; i++) {
        for(int j=0; j<ydim; j++) {
            qREAL ytmp = y[i*ydim+j];
            if(isnanq(ytmp) || isinfq(ytmp)) {
                self->nNAN++;
                if(self->nNAN > self->NANMax) {
                    NaNQ = true;
                    break;
                } else {
                    y[i*ydim+j] = 0;
                }
            }
        }
        if(self->ReIm == 1) {
            y[i*ydim+1] = 0;
        } else if(self->ReIm == 2) {
            y[i*ydim+0] = 0;
        }
    }
        
    return NaNQ ? 1 : 0;
}

void HCubature::DefaultPrintHooker(qREAL* result, qREAL* epsabs, long long int* nrun, void *fdata) {
    auto self = (HCubature*)fdata;
    if(*nrun == self->MaxPTS + 1979) return;

    if(self->Verbose>10 && self->RunMAX>0 && (*nrun-self->NEval) >= self->RunPTS ) {
        char r0[64], r1[64], e0[32], e1[32];
        quadmath_snprintf(r0, sizeof r0, "%.10QG", result[0]);
        quadmath_snprintf(r1, sizeof r1, "%.10QG", result[1]);
        quadmath_snprintf(e0, sizeof e0, "%.5QG", epsabs[0]);
        quadmath_snprintf(e1, sizeof e1, "%.5QG", epsabs[1]);
        cout << "     N: " << (*nrun) << ", ";
        if(self->ReIm==3 || self->ReIm==1) cout << "["<<r0 << ", " << e0 << "]";
        if(self->ReIm==3 || self->ReIm==2) cout << "+I*[" << r1 << ", " << e1 << "]";
        cout << endl;
    }
    self->NEval = *nrun;
    
    if((isnanq(result[0]) || isnanq(result[1]) || isnanq(epsabs[0]) || isnanq(epsabs[1])) || (isinfq(result[0]) || isinfq(result[1]) || isinfq(epsabs[0]) || isinfq(epsabs[1]))) {
         *nrun = self->MaxPTS + 1979;
         if(self->LastState>0) self->LastState = -1;
         if(self->Verbose>10 && self->RunMAX>0) cout << RED << "     Exit: NaN, N = " << self->NEval << RESET << endl;
         return;
    }
    
    if(self->RunMAX>0 && (epsabs[0] > 1E30*self->EpsAbs || epsabs[1] > 1E30*self->EpsAbs)) {
         *nrun = self->MaxPTS + 1979;
         if(self->LastState>0) self->LastState = -1;
         if(self->Verbose>10 && self->RunMAX>0) cout << RED << "     Exit: EpsAbs, N = " << self->NEval << RESET << endl;
         return;
    }
    
    if((self->LastState == 0) || (epsabs[0]<=2*self->LastAbsErr[0] && epsabs[1]<=2*self->LastAbsErr[1])) {
        self->LastResult[0] = result[0];
        self->LastResult[1] = result[1];
        self->LastAbsErr[0] = epsabs[0];
        self->LastAbsErr[1] = epsabs[1];
        self->LastState = 1;
        self->lastNRUN = *nrun;
        self->lastnNAN = self->nNAN;
    }
    
    bool rExit = (epsabs[0] < self->EpsAbs+1E-50Q) || (epsabs[0] < fabsq(result[0])*self->EpsRel+1E-50Q);
    bool iExit = (epsabs[1] < self->EpsAbs+1E-50Q) || (epsabs[1] < fabsq(result[1])*self->EpsRel+1E-50Q);
    if(rExit && iExit) {
        *nrun = self->MaxPTS + 1979;
        return;
    }

    auto pid = getpid();
    ostringstream fn;
    fn << pid << ".int.done";
    if(file_exists(fn.str().c_str())) {
        *nrun = self->MaxPTS + 1979;
        ostringstream cmd;
        cmd << "rm " << fn.str();
        system(cmd.str().c_str());
        if(self->Verbose>10) cout << "     Exit: " << fn.str() << endl;
    }
}

ex HCubature::Integrate() {
    assert(mpfr_buildopt_tls_p()>0);
    
    unsigned int xdim = XDim;
    unsigned int ydim = 2;
    qREAL result[ydim], estabs[ydim];

    qREAL xmin[xdim], xmax[xdim];
    for(int i=0; i<xdim; i++) {
        xmin[i] = 0.0Q;
        xmax[i] = 1.0Q;
    }
    LastState = 0;
    NEval = 0;
    nNAN = 0;
    
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
        NEval = lastNRUN;
        nNAN = lastnNAN;
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
    
    // Check nNAN / NEval
    if(nNAN * 1000 > NEval) {
        cout << RED << "NAN=" << nNAN << " v.s. RUN=" << NEval << RESET << endl;
    }
    
    return FResult;
}

}
