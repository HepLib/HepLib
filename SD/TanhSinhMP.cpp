/**
 * @file
 * @brief Numerical Integrator using TanhSinhMP
 */

#include <math.h>
#include <complex>
extern "C" {
#include <quadmath.h>
}
#include "mpreal.h"
#include <omp.h>
#include "SD.h"

/* error return codes */
#define SUCCESS 0
#define FAILURE 1

using namespace std;
namespace {
    typedef mpfr::mpreal mpREAL;
    typedef const mpREAL & mpREAL_t;
    typedef complex<mpREAL> mpCOMPLEX;
    typedef std::function<int(unsigned ydim, mpREAL *y, mpREAL *e, const mpREAL & x, void *fdata)> f1Type;
    typedef const f1Type & f1Type_t;
    typedef std::function<int(unsigned ydim, mpREAL *y, mpREAL *e, unsigned xdim, const mpREAL *x, void *fdata)> fnType;
    typedef const fnType & fnType_t;
    typedef void (*PrintHookerType) (mpREAL *, mpREAL *, size_t *, void *);
}
extern mpREAL mpPi;
extern mpREAL mpEuler;

namespace {
    int CPUCORES = 8;
    // int f from a to b, parallel version
    int TanhSinh1(mpREAL *oval, mpREAL *oerr, f1Type_t f, unsigned ydim, mpREAL_t epsrel, PrintHookerType PrintHooker, void *fdata) {
        mpREAL a=0, b=1; // integration domain (a,b)
        mpREAL c = (a+b)/2; // gamma
        mpREAL d = (b-a)/2; // sigma
        const mpREAL tol = 10*epsrel;
        mpREAL s[ydim], e_s[ydim], val[ydim], err[ydim], v[ydim], h = 2;
        if(f(ydim, s, e_s, c, fdata)) return FAILURE;
        size_t k = 0;
        long kmax = ((HepLib::SD::TanhSinhMP*)fdata)->LevelMAX;
        if(kmax<0) kmax = -kmax;
        while(k<=kmax) {
            mpREAL p[ydim], fp[ydim], fm[ydim], q, t, eh;
            for(int j=0; j<ydim; j++) p[j] = fp[j] = fm[j] = 0;
            h /= 2;
            eh = exp(h);
            t = eh;
            if (k > 0) eh *= eh;
            bool parallel = true;
            if(omp_in_parallel() && ( !omp_get_nested() || omp_get_active_level()>1 )) parallel = false;
            if(parallel) { // Parallel
                while(true) {
                    int total = CPUCORES/2;
                    mpREAL xs[total], ws[total], fps[total*ydim], e_fps[total*ydim], fms[total*ydim], e_fms[total*ydim];
                    int RC1[total], RC2[total];
                    for(int i=0; i<total; i++) {
                        mpREAL u = exp(1/t-t);
                        mpREAL r = 2*u/(1+u);
                        ws[i] = (t+1/t)*r/(1+u);
                        xs[i] = d*r;
                        t *= eh;
                        RC1[i] = RC2[i] = 0;
                    }
                    auto prec = mpfr::mpreal::get_default_prec();
                    auto rnd = mpfr::mpreal::get_default_rnd();
                    #pragma omp parallel for
                    for(int ii=0; ii<2*total; ii++) {
                        int i = ii/2, i2 = ii%2;
                        mpfr::mpreal::set_default_prec(prec);
                        mpfr::mpreal::set_default_rnd(rnd);
                        mpREAL_t x = xs[i];
                        if(i2==0) { if (a+x > a) RC1[i] = f(ydim,fps+ydim*i,e_fps+ydim*i,a+x,fdata); }
                        else if(i2==1) { if (b-x < b) RC2[i] = f(ydim,fms+ydim*i,e_fms+ydim*i,b-x,fdata); }
                        mpfr_free_cache();
                    }
                    bool ok;
                    for(int i=0; i<total; i++) {
                        if(RC1[i]!=0 || RC2[i]!=0) return FAILURE;
                        mpREAL_t x = xs[i];
                        ok = true;
                        for(int j=0; j<ydim; j++) {
                            if (a+x > a) if (isfinite(fps[i])) {
                                fp[j] = fps[i*ydim+j];
                                if(e_s[j] < e_fps[i*ydim+j]) e_s[j] = e_fps[i*ydim+j];
                            }
                            if (b-x < b) if (isfinite(fms[i])) {
                                fm[j] = fms[i*ydim+j];
                                if(e_s[j] < e_fms[i*ydim+j]) e_s[j] = e_fms[i*ydim+j];
                            }
                            q = ws[i]*(fp[j]+fm[j]);
                            p[j] += q;
                            ok = ok && (fabs(q)<=epsrel*fabs(p[j]));
                        }
                        if(ok) break;
                    }
                    if(ok) break;
                }
            } else { // Non-Parallel
                while(true) {
                    mpREAL u = exp(1/t-t);
                    mpREAL r = 2*u/(1+u);
                    mpREAL w = (t+1/t)*r/(1+u);
                    mpREAL x = d*r;
                    t *= eh;
                    mpREAL y[ydim], e_y[ydim];
                    if (a+x > a) {
                        if(f(ydim,y,e_y,a+x,fdata)) return FAILURE;
                        for(int j=0; j<ydim; j++) if (isfinite(y[j])) {
                            fp[j] = y[j];
                            if(e_s[j] < e_y[j]) e_s[j] = e_y[j];
                        }
                    }
                    if (b-x < b) {
                        if(f(ydim,y,e_y,b-x,fdata)) return FAILURE;
                        for(int j=0; j<ydim; j++) if (isfinite(y[j])) {
                            fm[j] = y[j];
                            if(e_s[j] < e_y[j]) e_s[j] = e_y[j];
                        }
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
                err[j] = fabs(d*h*v[j]) + (b-a)*e_s[j];
            }
            if(k > 10) {
                cout << RED << "Large Level: " << RESET << k << endl;
                for(int j=0; j<ydim; j++) {
                    cout << val[j] << " +- " << err[j].toString(3) << endl;
                }
                cout << endl;
            }
            size_t kk = k;
            if(PrintHooker) PrintHooker(val, err, &kk, fdata);
            if(kk!=k) {
                for(int j=0; j<ydim; j++) {
                    oval[j] = val[j];
                    oerr[j] = err[j];
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
            oval[j] = val[j];
            oerr[j] = err[j];
        }
        return SUCCESS;
    }
    
    int TanhSinhN(mpREAL *val, mpREAL *err, fnType_t f, unsigned xdim, unsigned ydim, mpREAL_t eps, PrintHookerType PrintHooker, void *fdata) {
        if(xdim==1) {
            auto f1 = [f](unsigned ydim, mpREAL *y, mpREAL *e, mpREAL_t x, void *fdata)->int {
                mpREAL xs[1];
                xs[0] = x;
                if(f(ydim,y,e,1,xs,fdata)) return FAILURE;
                return SUCCESS;
            };
            return TanhSinh1(val,err,f1,ydim,eps,PrintHooker,fdata);
        } else {
            auto f1 = [f,xdim,eps](unsigned ydim, mpREAL *y, mpREAL *e, mpREAL_t x, void *fdata)->int {
                auto f2 =[f,x](unsigned ydim, mpREAL *y1, mpREAL *e1, unsigned xdim1, const mpREAL *xs1, void *fdata)->int {
                    mpREAL xs[xdim1+1];
                    if(true) { // insert first
                        xs[0] = x;
                        for(int i=1; i<=xdim1; i++) xs[i] = xs1[i-1];
                    } else { // insert last
                        for(int i=0; i<xdim1; i++) xs[i] = xs1[i];
                        xs[xdim1] = x;
                    }
                    if(f(ydim,y1,e1,xdim1+1,xs,fdata)) return FAILURE;
                    return SUCCESS;
                };
                return TanhSinhN(y,e,f2,xdim-1,ydim,eps/xdim,NULL,fdata);
            };
            return TanhSinh1(val,err,f1,ydim,eps,PrintHooker,fdata);
        }
    }
}

namespace HepLib::SD {

/*-----------------------------------------------------*/
// TanhSinhMP Classes
/*-----------------------------------------------------*/

TanhSinhMP::TanhSinhMP(size_t kmax) { LevelMAX = kmax; }

ex TanhSinhMP::mp2ex(mpREAL_t num) {
    ostringstream oss;
    oss.precision(MPDigits);
    oss << num;
    string sn = oss.str();
    numeric ret(sn.c_str());
    return ret;
}

int TanhSinhMP::Wrapper(unsigned ydim, mpREAL *y, mpREAL *e, unsigned xdim, const mpREAL *x, void *fdata) {
    for(int j=0; j<ydim; j++) e[j] = 0; // assume error can be ignored
    auto self = (TanhSinhMP*)fdata;
    self->IntegrandMP(xdim, x, ydim, y, self->mpParameter, self->mpLambda);
    bool ok = true;
    for(int j=0; j<ydim; j++) {
        if(!isfinite(y[j])) { ok = false; break; }
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
    if(*nrun == self->LevelMAX + 1979) return;
    if(self->RunTime>0) {
        auto cur_timer = time(NULL);
        auto used_time = difftime(cur_timer,self->StartTimer);
        if(used_time>self->RunTime) {
            self->NEval = *nrun;
            *nrun = self->LevelMAX + 1979;
            if(Verbose>10) cout << WarnColor << "     Exit with Run out of Time: " << used_time << RESET << endl;
            return;
        }
    }
    if(Verbose>10 && self->LevelMAX>0) {
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
        *nrun = self->LevelMAX + 1979;
        return;
    }

    auto pid = getpid();
    ostringstream fn;
    fn << pid << ".int.done";
    if(file_exists(fn.str().c_str())) {
        self->NEval = *nrun;
        *nrun = self->LevelMAX + 1979;
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
    if(!isfinite(result[0]) || !isfinite(result[1])) FResult += NaN;
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
