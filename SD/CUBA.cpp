/**
 * @file
 * @brief Numerical Integrator using CUBA
 */
 
#include "SD.h"
#include <math.h>
#include <complex>
extern "C" {
#include <quadmath.h>
}
#include "mpreal.h"

using namespace std;
typedef __float128 qREAL;
typedef __complex128 qCOMPLEX;
typedef long double dREAL;
typedef complex<dREAL> dCOMPLEX;
typedef mpfr::mpreal mpREAL;
typedef complex<mpREAL> mpCOMPLEX;

extern const qCOMPLEX qiEpsilon;
extern mpREAL mpPi;
extern mpREAL mpEuler;
extern mpCOMPLEX mpiEpsilon;

namespace HepLib::SD {

/*-----------------------------------------------------*/
// CUBA Classes
/*-----------------------------------------------------*/
extern "C" {
#include "cubaq.h"
}

int CUBA::Wrapper(const int *pxdim, const qREAL *x, const int *pydim, qREAL *y, void *fdata) {
    auto self = (CUBA*)fdata;
    int xdim = *pxdim, ydim = *pydim;
    bool NaNQ = false;
    //#pragma omp parallel for num_threads(omp_get_num_procs()-1) schedule(dynamic, 1)
    int npts = 1;
    for(int i=0; i<npts; i++) {
        int iDQMP = self->inDQMP(x+i*xdim);
        if((self->IntegrandMP!=NULL) && (self->DQMP>2 || iDQMP>2)) {
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
        if(!ok) {
            qREAL xx[xdim];
            for(int ii=0; ii<xdim; ii++) xx[ii] = x[i*xdim+ii] * 0.995Q;
            self->IntegrandMP(xdim, xx, ydim, y+i*ydim, self->Parameter, self->Lambda);
        }
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

ex CUBA::Integrate() {
    mpfr_free_cache();
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(MPDigits));
    mpPi = mpfr::const_pi();
    mpiEpsilon = complex<mpREAL>(0,cimagq(qiEpsilon));
    mpEuler = mpfr::const_euler();
    if(mpfr_buildopt_tls_p()<=0) {
        throw Error("Integrate: mpfr_buildopt_tls_p()<=0.");
    }
        
    unsigned int xdim = XDim;
    unsigned int ydim = 2;
    qREAL result[ydim], estabs[ydim];
    
    MaxPTS = RunPTS * RunMAX;
    if(MaxPTS<0) MaxPTS = -MaxPTS;
    
    int ndim = xdim;
    long long int NVEC = 1;
    
    int nregions, fail;
    long long neval;
    qREAL prob[ydim];
    
    switch(Method) {
        case VEGAS:
            llVegas(ndim, ydim, Wrapper, this, NVEC, EpsRel, EpsAbs, VERBOSE,
                VEGAS_SEED, RunPTS, MaxPTS, VEGAS_NSTART, VEGAS_NINCREASE, VEGAS_NBATCH,
                0, NULL, NULL,
                &neval, &fail, result, estabs, prob
            );
            break;
        case CUHRE:
            if(ndim==1) ndim = 2;
            llCuhre(ndim, ydim, Wrapper, this, NVEC, EpsRel, EpsAbs, VERBOSE | 4, RunPTS, MaxPTS,
                CUHRE_KEY, NULL, NULL,
                &nregions, &neval, &fail, result, estabs, prob
            );
            break;
    }
    NEval = neval;
    
    ex FResult = 0;
    if(isnanq(result[0]) || isnanq(result[1])) FResult += NaN;
    else {
        try{
            FResult += VE(q2ex(result[0]), q2ex(estabs[0]));
            FResult += VE(q2ex(result[1]), q2ex(estabs[1])) * I;
        } catch(...) {
            FResult += NaN;
        }
    }
    
    return FResult;
}

}
