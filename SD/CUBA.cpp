#include "SD.h"

namespace HepLib {

/*-----------------------------------------------------*/
// CUBA Classes
/*-----------------------------------------------------*/
extern "C" {
#include "cubaq.h"
}

bool CUBA::useQ(unsigned xdim, qREAL const *x) {
    if(xdim<2) return true;
    for(int i=0; i<xdim; i++) {
        if(x[i] < 5.0E-5Q) return true;
    }
    return false;
}

int CUBA::Wrapper(const int *pxdim, const qREAL *x, const int *pydim, qREAL *y, void *fdata) {
    auto self = (CUBA*)fdata;
    int xdim = *pxdim, ydim = *pydim;
    bool NaNQ = false;
    //#pragma omp parallel for num_threads(omp_get_num_procs()) schedule(dynamic, 1)
    int npts = 1;
    for(int i=0; i<npts; i++) {
        if(self->UseQ || useQ(xdim, x+i*xdim)) {
            self->IntegrandQ(xdim, x+i*xdim, ydim, y+i*ydim, self->Parameter, self->Lambda);
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

    }
    return NaNQ ? 1 : 0;
}

ex CUBA::Integrate(unsigned int xdim, SD_Type fp, SD_Type fpQ, const qREAL* pl, const qREAL* la) {
    Integrand = fp;
    IntegrandQ = fpQ;
    Parameter = pl;
    Lambda = la;
    
    int ydim = 2;
    qREAL result[ydim], estabs[ydim];
    
    MaxPTS = RunPTS * RunMAX;
    if(MaxPTS<0) MaxPTS = -MaxPTS;
    
    int ndim = xdim;
    long long int nvec = 1;
    int VERBOSE = Verbose;
    
    int nregions, fail;
    long long neval;
    qREAL prob[ydim];
    
    switch(Method) {
        case VEGAS:
            llVegas(ndim, ydim, Wrapper, this, nvec, EpsRel, EpsAbs, VERBOSE,
                VEGAS_SEED, RunPTS, MaxPTS, VEGAS_NSTART, VEGAS_NINCREASE, VEGAS_NBATCH,
                0, NULL, NULL,
                &neval, &fail, result, estabs, prob
            );
            break;
        case CUHRE:
            if(ndim==1) ndim = 2;
            llCuhre(ndim, ydim, Wrapper, this, nvec, EpsRel, EpsAbs, VERBOSE | 4, RunPTS, MaxPTS,
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
