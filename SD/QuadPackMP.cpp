/**
 * @file
 * @brief Numerical Integrator using QuadPackMP
 */
 
#include "SD.h"
#include <math.h>
#include <complex>
extern "C" {
#include <quadmath.h>
}
#include "mpreal.h"
#include "Lib3_QuadPack.h"

/* error return codes */
#define SUCCESS 0
#define FAILURE 1

using namespace std;
namespace {
    typedef mpfr::mpreal mpREAL;
    typedef const mpREAL & mpREAL_t;
    typedef complex<mpREAL> mpCOMPLEX;
    typedef std::function<int(mpREAL &y, mpREAL &e, const mpREAL & x, void *fdata)> f1Type;
    typedef std::function<int(mpREAL &y, mpREAL &e, unsigned xdim, const mpREAL *x, void *fdata)> fnType;
    typedef void (*PrintHookerType) (mpREAL *, mpREAL *, size_t *, void *);
}

extern mpREAL mpPi;
extern mpREAL mpEuler;
extern mpCOMPLEX mpiEpsilon;

namespace {
    int CPUCORES = 8;
    size_t QAG_n = 10000;
    size_t QAG_m = 10; // sets (2m+1)-point Gauss-Kronrod
    
    int QuadPack1(mpREAL &oval, mpREAL &oerr, f1Type f, mpREAL_t epsabs, PrintHookerType PrintHooker, void *fdata) {
        Workspace Work(QAG_n, QAG_m);
        mpREAL result, abserr;
        Function F(f,fdata);
        try {
            return Work.qag(F, 0, 1, epsabs, 0, oval, oerr, PrintHooker);
        } catch (const char* reason) {
            cout << reason << endl;
            throw reason;
        }
        return SUCCESS;
    }
    
    mpREAL QuadPackN(mpREAL * oerr, const int & xdim, const fnType & f, const mpREAL & eps, PrintHookerType PrintHooker, void * fdata) {
        if(xdim==1) {
            auto f1 = [f,eps](const mpREAL & x, void * fdata)->mpREAL {
                vector<mpREAL> xs;
                xs.push_back(x);
                return f(xs,fdata);
            };
            return QuadPack1(oerr, f1, eps, PrintHooker, fdata);
        }
        
        auto f1 = [f,eps,xdim](const mpREAL & x, void * fdata)->mpREAL {
            auto f2 =[f,x,xdim](const vector<mpREAL> & xs1, void * fdata1)->mpREAL {
                vector<mpREAL> xs(xdim);
                if(true) { // insert first
                    xs[0] = x;
                    for(int i=1; i<xdim; i++) xs[i] = xs1[i-1];
                } else { // insert last
                    for(int i=0; i<xdim-1; i++) xs[i] = xs1[i-1];
                    xs[xdim-1] = x;
                }
                return f(xs,fdata1);
            };
            return QuadPackN(NULL, xdim-1, f2, eps/xdim, NULL, fdata);
        };
        return QuadPack1(oerr, f1, eps, PrintHooker, fdata);
    }
    
}

namespace HepLib::SD {

/*-----------------------------------------------------*/
// QuadPackMP Classes
/*-----------------------------------------------------*/

ex QuadPackMP::mp2ex(const mpREAL & num) {
    ostringstream oss;
    oss.precision(MPDigits);
    oss << num;
    string sn = oss.str();
    numeric ret(sn.c_str());
    return ret;
}

mpREAL QuadPackMP::Wrapper(const vector<mpREAL> & xs, void * fdata) {

    auto self = (QuadPackMP*)fdata;
    int xdim = xs.size(), ydim = 2;
    mpREAL x[xdim], y[ydim];
    for(int i=0; i<xdim; i++) x[i] = xs[i];
    self->IntegrandMP(xdim, x, ydim, y, self->mpParameter, self->mpLambda);

    // Final Check NaN/Inf
    bool ok = true;
    for(int j=0; j<ydim; j++) {
        mpREAL ytmp = y[j];
        if(isnan(ytmp) || isinf(ytmp)) { ok = false; break; }
    }
    if(!ok) {
        mpfr::mpreal::set_default_prec(mpfr::digits2bits(self->MPDigits*100));
        self->IntegrandMP(xdim, x, ydim, y, self->mpParameter, self->mpLambda);
        mpfr::mpreal::set_default_prec(mpfr::digits2bits(self->MPDigits));
    }
    
    // final check
    for(int j=0; j<ydim; j++) {
        mpREAL ytmp = y[j];
        if(isnan(ytmp) || isinf(ytmp)) {
            #pragma omp atomic
            self->nNAN++;
            if(self->nNAN > self->NANMax) break;
            else y[j] = 0;
        }
    }

    return y[self->Index];
}

void QuadPackMP::DefaultPrintHooker(mpREAL* result, mpREAL* epsabs, size_t * nrun, void *fdata) {
    auto self = (QuadPackMP*)fdata;
    if(*nrun == self->RunMAX + 1979) return;
    if(self->RunTime>0) {
        auto cur_timer = time(NULL);
        auto used_time = difftime(cur_timer,self->StartTimer);
        if(used_time>self->RunTime) {
            self->NEval = *nrun;
            *nrun = self->RunMAX + 1979;
            if(Verbose>10) cout << WarnColor << "     Exit with Run out of Time: " << used_time << RESET << endl;
            return;
        }
    }
    if(Verbose>10) {
        auto r0 = result[0];
        auto r1 = result[1];
        auto e0 = epsabs[0].toString(3);
        auto e1 = epsabs[1].toString(3);
        cout << "     L: " << (*nrun) << ", ";
        if(self->ReIm==3 || self->ReIm==1) cout << "[" << r0 << ", " << e0 << "]";
        if(self->ReIm==3 || self->ReIm==2) cout << "+I*[" << r1 << ", " << e1 << "]";
        cout << endl;
    }
    self->NEval = *nrun;
    
    if((isnan(result[0]) || isnan(result[1]) || isnan(epsabs[0]) || isnan(epsabs[1])) || (isinf(result[0]) || isinf(result[1]) || isinf(epsabs[0]) || isinf(epsabs[1]))) {
        self->NEval = *nrun;
        *nrun = self->RunMAX + 1979;
        if(self->LastState>0) self->LastState = -1;
        if(Verbose>10) cout << ErrColor << "     Exit with NaN, LastN=" << self->lastNRUN << RESET << endl;
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

    bool rExit = (epsabs[0] < self->EpsAbs+1E-50Q) || (epsabs[0] < fabs(result[0])*self->EpsRel+1E-50Q);
    bool iExit = (epsabs[1] < self->EpsAbs+1E-50Q) || (epsabs[1] < fabs(result[1])*self->EpsRel+1E-50Q);
    if(rExit && iExit) {
        *nrun = self->RunMAX + 1979;
        return;
    }

    auto pid = getpid();
    ostringstream fn;
    fn << pid << ".int.done";
    if(file_exists(fn.str().c_str())) {
        self->NEval = *nrun;
        *nrun = self->RunMAX + 1979;
        ostringstream cmd;
        cmd << "rm " << fn.str();
        system(cmd.str().c_str());
        if(Verbose>10) cout << "     Exit: " << fn.str() << endl;
    }
}

ex QuadPackMP::Integrate() {
    CPUCORES = omp_get_num_procs();
    if(mpfr_buildopt_tls_p()<=0) throw Error("Integrate: mpfr_buildopt_tls_p()<=0.");
    mpfr_free_cache();
    mpfr::mpreal::set_default_prec(mpfr::digits2bits(MPDigits));
    mpPi = mpfr::const_pi();
    mpEuler = mpfr::const_euler();
    mpiEpsilon = complex<mpREAL>(0,mpiEpsilon.imag());
    
    unsigned int xdim = XDim;
    unsigned int ydim = 2;
    mpREAL result[ydim], estabs[ydim];

    LastState = 0;
    NEval = 0;
    nNAN = 0;
    
    QAG_n = nQAG;
    QAG_m = mQAG;
    StartTimer = time(NULL);
    Index = 0;
    result[Index] = QuadPackN(estabs+Index, xdim, Wrapper, EpsAbs, PrintHooker, this);
    StartTimer = time(NULL);
    Index = 1;
    result[Index] = QuadPackN(estabs+Index, xdim, Wrapper, EpsAbs, PrintHooker, this);
    
    
    int nok;// = QuadPackN(result, estabs, Wrapper, xdim, ydim, EpsAbs, PrintHooker, this);

    if(nok) {
        mpREAL abs_res = sqrt(result[0]*result[0]+result[1]*result[1]);
        mpREAL abs_est = sqrt(estabs[0]*estabs[0]+estabs[1]*estabs[1]);
        mpREAL mpfr_eps = 10*mpfr::machine_epsilon();
        if( (abs_res < mpfr_eps) && (abs_est < mpfr_eps) ) {
            cout << ErrColor << "QuadPackMP Failed with 0 result returned!" << RESET << endl;
            return NaN;
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
    if(isnan(result[0]) || isnan(result[1])) FResult += NaN;
    else {
        try{
            FResult += VE(mp2ex(result[0]), mp2ex(estabs[0]));
            FResult += VE(mp2ex(result[1]), mp2ex(estabs[1])) * I;
        } catch(...) {
            FResult += NaN;
        }
    }

    return FResult;
}

}
