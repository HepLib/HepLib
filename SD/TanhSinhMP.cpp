/**
 * @file
 * @brief Numerical Integrator using TanhSinhMP
 */
 
#include "SD.h"
#include <math.h>
#include <complex>
extern "C" {
#include <quadmath.h>
}
#include "mpreal.h"

/* error return codes */
#define SUCCESS 0
#define FAILURE 1

using namespace std;
typedef mpfr::mpreal mpREAL;
typedef complex<mpREAL> mpCOMPLEX;

extern mpREAL mpPi;
extern mpREAL mpEuler;
extern mpCOMPLEX mpiEpsilon;

using namespace std;
typedef std::function<int(unsigned ydim, mpREAL * y, const mpREAL & x, void * fdata)> f1Type;
typedef std::function<int(unsigned ydim, mpREAL * y, unsigned xdim, const mpREAL * x, void * fdata)> fnType;
typedef void (*PrintHookerType) (mpREAL*, mpREAL*, size_t *, void *);

namespace {
    int CPUCORES = 8;
    // int f from a to b, parallel version
    int TanhSinh1(mpREAL *oval, mpREAL *oerr, f1Type f, unsigned ydim, const mpREAL & sign_epsrel, const PrintHookerType & PrintHooker=NULL, void* fdata=NULL) {
        mpREAL epsrel = sign_epsrel;
        bool parallel = false;
        if(epsrel<0) {
            epsrel = 0-epsrel;
            parallel = true;
        }
        mpREAL a=0, b=1; // integration domain (a,b)
        mpREAL c = (a+b)/2; // gamma
        mpREAL d = (b-a)/2; // sigma
        const mpREAL tol = 10*epsrel;
        mpREAL s[ydim], val[ydim], err[ydim], v[ydim], h = 2;
        if(f(ydim, s, c, fdata)) return FAILURE;
        size_t k = 0;
        size_t kmax = ((HepLib::SD::TanhSinhMP*)fdata)->RunMAX;
        while(k<=kmax) {
            mpREAL p[ydim], fp[ydim], fm[ydim], q, t, eh;
            for(int j=0; j<ydim; j++) p[j] = fp[j] = fm[j] = 0;
            h /= 2;
            eh = exp(h);
            t = eh;
            if (k > 0) eh *= eh;
            if(parallel) {
                while(true) {
                    int total = CPUCORES/2; // TODO: tune this parameter
                    mpREAL xs[total], ws[total], fps[total*ydim], fms[total*ydim];
                    int RC[total];
                    for(int i=0; i<total; i++) {
                        mpREAL u = exp(1/t-t);
                        mpREAL r = 2*u/(1+u);
                        ws[i] = (t+1/t)*r/(1+u);
                        xs[i] = d*r;
                        t *= eh;
                        RC[i] = 0;
                    }
                    auto prec = mpfr::mpreal::get_default_prec();
                    auto rnd = mpfr::mpreal::get_default_rnd();
                    #pragma omp parallel for schedule(dynamic,1)
                    for(int i=0; i<total; i++) {
                        mpfr::mpreal::set_default_prec(prec);
                        mpfr::mpreal::set_default_rnd(rnd);
                        const mpREAL & x = xs[i];
                        if (a+x > a) RC[i]+=f(ydim,fps+ydim*i,a+x,fdata);
                        if (b-x < b) RC[i]+f(ydim,fms+ydim*i,b-x,fdata);
                        mpfr_free_cache();
                    }
                    bool ok;
                    for(int i=0; i<total; i++) {
                        if(RC[i]!=0) return FAILURE;
                        const mpREAL & x = xs[i];
                        ok = true;
                        for(int j=0; j<ydim; j++) {
                            if (a+x > a) if (isfinite(fps[i])) fp[j] = fps[i*ydim+j];
                            if (b-x < b) if (isfinite(fms[i])) fm[j] = fms[i*ydim+j];
                            q = ws[i]*(fp[j]+fm[j]);
                            p[j] += q;
                            ok = ok && (fabs(q)<=epsrel*fabs(p[j]));
                        }
                        if(ok) break;
                    }
                    if(ok) break;
                }
            } else {
                while(true) {
                    mpREAL u = exp(1/t-t);
                    mpREAL r = 2*u/(1+u);
                    mpREAL w = (t+1/t)*r/(1+u);
                    mpREAL x = d*r;
                    t *= eh;
                    mpREAL y[ydim];
                    if (a+x > a) {
                        if(f(ydim,y,a+x,fdata)) return FAILURE;
                        for(int j=0; j<ydim; j++) { if (isfinite(y[j])) fp[j] = y[j]; }
                    }
                    if (b-x < b) {
                        if(f(ydim,y,b-x,fdata)) return FAILURE;
                        for(int j=0; j<ydim; j++) { if (isfinite(y[j])) fm[j] = y[j]; }
                    }
                    bool ok = true;
                    for(int j=0; j<ydim; j++) {
                        q = w*(fp[j]+fm[j]);
                        p[j] += q;
                        ok = ok && (fabs(q)<=epsrel*fabs(p[j]));
                    }
                    if(ok) break;
                }
            }
            for(int j=0; j<ydim; j++) {
                v[j] = s[j] - p[j];
                s[j] += p[j];
                val[j] = d*h*s[j];
                err[j] = fabs(d*h*v[j]);
                if(fabs(val[j]) > 1E10) cout << RED << "Large fValue: " << RESET << val[j] << endl;
            }
            if(k > 10) cout << RED << "Large Level: " << RESET << k << endl;
            size_t kk = k;
            if(PrintHooker) PrintHooker(val, err, &kk, fdata);
            if(kk!=k) {
                for(int j=0; j<ydim; j++) {
                    if(oval!=NULL) oval[j] = val[j];
                    if(oerr!=NULL) oerr[j] = err[j];
                }
                return SUCCESS;
            }
            bool ok = true;
            for(int j=0; j<ydim; j++) {
                ok = ok && (fabs(v[j])<=tol*fabs(s[j]));
                if(!ok) break;
            }
            if(ok) break;
            ++k;
        }
        for(int j=0; j<ydim; j++) {
            if(oval!=NULL) oval[j] = val[j];
            if(oerr!=NULL) oerr[j] = err[j];
        }
        return SUCCESS;
    }
        
    int TanhSinhN(mpREAL *val, mpREAL *err, fnType f, unsigned xdim, unsigned ydim, const mpREAL & eps, const PrintHookerType & PrintHooker=NULL, void* fdata=NULL) {
        if(xdim==1) {
            auto f1 = [f](unsigned ydim, mpREAL * y, const mpREAL & x, void * fdata)->int {
                mpREAL xs[1];
                xs[0] = x;
                if(f(ydim,y,1,xs,fdata)) return FAILURE;
                return SUCCESS;
            };
            return TanhSinh1(val,err,f1,ydim,-eps,PrintHooker,fdata); // eps<0 for parallel mode
        }
        
        auto f1 = [f,xdim,eps](unsigned ydim, mpREAL * y, const mpREAL & x, void * fdata)->int {
            auto f2 =[f,x](unsigned ydim, mpREAL * y1, unsigned xdim1, const mpREAL * xs1, void * fdata)->int {
                mpREAL xs[xdim1+1];
                for(int i=0; i<xdim1; i++) xs[i] = xs1[i];
                xs[xdim1] = x;
                if(f(ydim,y1,xdim1+1,xs,fdata)) return FAILURE;
                return SUCCESS;
            };
            return TanhSinhN(y,NULL,f2,xdim-1,ydim,eps/xdim,NULL,fdata);
        };
        return TanhSinh1(val,err,f1,ydim,eps,PrintHooker,fdata);
    }
}

namespace HepLib::SD {

/*-----------------------------------------------------*/
// TanhSinhMP Classes
/*-----------------------------------------------------*/

TanhSinhMP::TanhSinhMP(size_t kmax) { RunMAX = kmax; }

ex TanhSinhMP::mp2ex(const mpREAL & num) {
    ostringstream oss;
    oss.precision(MPDigits);
    oss << num;
    string sn = oss.str();
    numeric ret(sn.c_str());
    return ret;
}

int TanhSinhMP::Wrapper(unsigned ydim, mpREAL * y, unsigned xdim, const mpREAL * x, void * fdata) {

    auto self = (TanhSinhMP*)fdata;
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

    return SUCCESS;
}

void TanhSinhMP::DefaultPrintHooker(mpREAL* result, mpREAL* epsabs, size_t * nrun, void *fdata) {
    auto self = (TanhSinhMP*)fdata;
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

ex TanhSinhMP::Integrate() {
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
    NEval = 0;
    mpREAL epsrel = 1E-2;
    int nok = TanhSinhN(result, estabs, Wrapper, xdim, ydim, epsrel, NULL, this);
    epsrel = min(EpsAbs/(fabs(result[0])+EpsAbs), EpsAbs/(fabs(result[1])+EpsAbs));
    if(epsrel < 1E-2) {
        StartTimer = time(NULL);
        nok = TanhSinhN(result, estabs, Wrapper, xdim, ydim, epsrel, PrintHooker, this);
    }

    if(nok) {
        mpREAL abs_res = sqrt(result[0]*result[0]+result[1]*result[1]);
        mpREAL abs_est = sqrt(estabs[0]*estabs[0]+estabs[1]*estabs[1]);
        mpREAL mpfr_eps = 10*mpfr::machine_epsilon();
        if( (abs_res < mpfr_eps) && (abs_est < mpfr_eps) ) {
            cout << ErrColor << "TanhSinhMP Failed with 0 result returned!" << RESET << endl;
            return NaN;
        }
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
