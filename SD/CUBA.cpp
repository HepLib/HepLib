#include "SD.h"

namespace HepLib {

/*-----------------------------------------------------*/
// CUBA Classes
/*-----------------------------------------------------*/
extern "C" {
#include "cubaq.h"
}

int CUBA::inDQMP(unsigned xdim, qREAL const *x) {
    return 3;
    if(xdim<2) return 3;
    for(int i=0; i<xdim; i++) {
        if(x[i] < 5.0E-5Q) return 3;
    }
    return 1;
}

int CUBA::Wrapper(const int *pxdim, const qREAL *x, const int *pydim, qREAL *y, void *fdata) {
    auto self = (CUBA*)fdata;
    int xdim = *pxdim, ydim = *pydim;
    bool NaNQ = false;
    //#pragma omp parallel for num_threads(omp_get_num_procs()) schedule(dynamic, 1)
    int npts = 1;
    for(int i=0; i<npts; i++) {
        int iDQMP = inDQMP(xdim, x+i*xdim);
        if(self->DQMP>2 || iDQMP>2) {
            self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
        } else if(self->DQMP>1 || iDQMP>1) {
            self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            bool ok = true;
            for(int j=0; j<ydim; j++) {
                qREAL ytmp = y[i*ydim+j];
                if(isnanq(ytmp)) {
                    ok = false;
                    break;
                }
            }
            if(!ok) {
                self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            }
        } else {
            self->Integrand(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            bool ok = true;
            
            for(int j=0; j<ydim; j++) {
                qREAL ytmp = y[i*ydim+j];
                if(isnanq(ytmp)) {
                    ok = false;
                    break;
                }
            }
            if(!ok) {
                self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            }
            
            for(int j=0; j<ydim; j++) {
                qREAL ytmp = y[i*ydim+j];
                if(isnanq(ytmp)) {
                    ok = false;
                    break;
                }
            }
            if(!ok) {
                self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
            }
        }
    }
    return NaNQ ? 1 : 0;
}

ex CUBA::Integrate(unsigned int xdim, SD_Type fp, SD_Type fpQ, SD_Type fpMP, const qREAL* pl, const qREAL* la) {
    Integrand = fp;
    IntegrandQ = fpQ;
    IntegrandMP = fpMP;
    Parameter = pl;
    Lambda = la;
    
    int ydim = 2;
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
