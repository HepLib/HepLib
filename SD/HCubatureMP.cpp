/**
 * @file
 * @brief Numerical Integrator using HCubature
 */
 
#include <math.h>
#include <complex>
extern "C" {
#include <quadmath.h>
}
#include "mpreal.h"
#include "SD.h"

using namespace std;
typedef mpfr::mpreal mpREAL;
typedef complex<mpREAL> mpCOMPLEX;

extern mpREAL mpPi;
extern mpREAL mpEuler;
extern mpCOMPLEX mpiEpsilon;

#include "Lib3_HCubatureMP.h"
namespace HepLib::SD {

    ex HCubatureMP::mp2ex(const mpREAL & num) {
        ostringstream oss;
        oss.precision(MPDigits);
        oss << num;
        string sn = oss.str();
        numeric ret(sn.c_str());
        return ret;
    }

    int HCubatureMP::Wrapper(unsigned int xdim, size_t npts, const mpREAL *x, void *fdata, unsigned int ydim, mpREAL *y) {
        auto self = (HCubatureMP*)fdata;
        bool NaNQ = false;
        
        unsigned int nthreads = self->Threads>0 ? self->Threads : omp_get_num_procs();
        #pragma omp parallel for num_threads(nthreads) schedule(dynamic, 1)
        for(int i=0; i<npts; i++) {
            mpfr_free_cache();
            auto pbit = mpfr::digits2bits(self->MPDigits);
            if(mpfr::mpreal::get_default_prec()!=pbit) mpfr::mpreal::set_default_prec(pbit);
            self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->mpParameter, self->mpLambda);
            // Final Check NaN/Inf
            bool ok = true;
            for(int j=0; j<ydim; j++) {
                mpREAL ytmp = y[i*ydim+j];
                if(isnan(ytmp) || isinf(ytmp)) { ok = false; break; }
            }
            if(!ok && (self->IntegrandMP!=NULL)) {
                mpfr_free_cache();
                auto pbit = mpfr::digits2bits(self->MPDigits*10);
                if(mpfr::mpreal::get_default_prec()!=pbit) mpfr::mpreal::set_default_prec(pbit);
                self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->mpParameter, self->mpLambda);
                pbit = mpfr::digits2bits(self->MPDigits);
                if(mpfr::mpreal::get_default_prec()!=pbit) mpfr::mpreal::set_default_prec(pbit);
            }
            
            // final check
            for(int j=0; j<ydim; j++) {
                mpREAL ytmp = y[i*ydim+j];
                if(isnan(ytmp) || isinf(ytmp)) {
                    #pragma omp atomic
                    self->nNAN++;
                    if(self->nNAN > self->NANMax) { NaNQ = true; break; }
                    else y[i*ydim+j] = 0;
                }
            }
            if(self->YDim==2) {
                if(self->ReIm == 1) y[i*ydim+1] = 0;
                else if(self->ReIm == 2) y[i*ydim+0] = 0;
            }
            mpfr_free_cache();
        }

        return NaNQ ? 1 : 0;
    }

    void HCubatureMP::DefaultPrintHooker(mpREAL* result, mpREAL* epsabs, size_t * nrun, void *fdata) {
        auto self = (HCubatureMP*)fdata;
        if(*nrun == self->MaxPTS + 1979) return;
        if(self->RunTime>0) {
            auto cur_timer = time(NULL);
            auto used_time = difftime(cur_timer,self->StartTimer);
            if(used_time>self->RunTime) {
                self->NEval = *nrun;
                *nrun = self->MaxPTS + 1979;
                if(Verbose>10) cout << WarnColor << "     Exit with Run out of Time: " << used_time << RESET << endl;
                return;
            }
        }
        if(Verbose>10 && (*nrun-self->NEval) >= self->RunPTS) {
            if(self->YDim==2) {
                auto r0 = result[0];
                auto r1 = result[1];
                auto e0 = epsabs[0].toString(3);
                auto e1 = epsabs[1].toString(3);
                cout << "     N: " << (*nrun) << ", ";
                if(self->ReIm==3 || self->ReIm==1) cout << "[" << r0 << ", " << e0 << "]";
                if(self->ReIm==3 || self->ReIm==2) cout << "+I*[" << r1 << ", " << e1 << "]";
                cout << endl;
            } else {
                for(int j=0; j<self->YDim; j++) {
                    auto r = result[j];
                    auto e = epsabs[j].toString(3);
                    if(j==0) cout << "     N: " << (*nrun) << ", ";
                    else {
                        cout << "          ";
                        int ct = to_string(*nrun).length();
                        for(int k=0; k<ct; k++) cout << " ";
                    }
                    cout << "[" << r << ", " << e << "]";
                    cout << endl;
                }
            }
        }
        if((*nrun-self->NEval) >= self->RunPTS) self->NEval = *nrun;
        
        bool any_NAN_Inf = false;
        for(int j=0; j<self->YDim; j++) {
            if(isnan(result[j]) || isnan(epsabs[j]) || isinf(result[j]) || isinf(epsabs[j])) {
                any_NAN_Inf = true;
                break;
            }
        }
        if(any_NAN_Inf) {
            self->NEval = *nrun;
            *nrun = self->MaxPTS + 1979;
            if(self->LastState>0) self->LastState = -1;
            if(Verbose>10) cout << ErrColor << "     Exit with NaN, LastN=" << self->lastNRUN << RESET << endl;
            return;
        }
        
        bool any_Big = false;
        for(int j=0; j<self->YDim; j++) {
            if(epsabs[j] > 1E30*self->EpsAbs) {
                any_Big = true;
                break;
            }
        }
        if(any_Big) {
            self->NEval = *nrun;
            *nrun = self->MaxPTS + 1979;
            if(self->LastState>0) self->LastState = -1;
            if(Verbose>10) cout << WarnColor << "     Exit with EpsAbs, LastN=" << self->lastNRUN << RESET << endl;
            return;
        }
        
        bool all_Done = true;
        for(int j=0; j<self->YDim; j++) {
            if(epsabs[j] > 2*self->LastAbsErr[j] ) {
                all_Done = false;
                break;
            }
        }
        if((self->LastState == 0) || all_Done) {
            for(int j=0; j<self->YDim; j++) {
                self->LastResult[j] = result[j];
                self->LastAbsErr[j] = epsabs[j];
            }
            self->LastState = 1;
            self->lastNRUN = *nrun;
            self->lastnNAN = self->nNAN;
        }
        
        bool bExit = true;
        for(int j=0; j<self->YDim; j++) {
            if((epsabs[j] > self->EpsAbs+1E-50Q) && (epsabs[j] > fabs(result[j])*self->EpsRel+1E-50Q)) {
                bExit = false;
                break;
            }
        }
        if(bExit && (*nrun)>self->MinPTS) {
            self->NEval = *nrun;
            *nrun = self->MaxPTS + 1979;
            return;
        }

        auto pid = getpid();
        ostringstream fn;
        fn << pid << ".int.done";
        if(file_exists(fn.str().c_str())) {
            self->NEval = *nrun;
            *nrun = self->MaxPTS + 1979;
            ostringstream cmd;
            cmd << "rm " << fn.str();
            auto rc = system(cmd.str().c_str());
            if(Verbose>10) cout << "     Exit: " << fn.str() << endl;
        }
    }

    ex HCubatureMP::Integrate(size_t tn) {
        if(LastResult.size()!=YDim) LastResult.resize(YDim);
        if(LastAbsErr.size()!=YDim) LastAbsErr.resize(YDim);
        if(mpfr_buildopt_tls_p()<=0) throw Error("Integrate: mpfr_buildopt_tls_p()<=0.");
        mpfr_free_cache();
        auto pbit = mpfr::digits2bits(MPDigits);
        if(mpfr::mpreal::get_default_prec()!=pbit) mpfr::mpreal::set_default_prec(pbit);
        mpPi = mpfr::const_pi();
        mpEuler = mpfr::const_euler();
        mpiEpsilon = complex<mpREAL>(0,mpfr::machine_epsilon()*100);
        
        unsigned int xdim = XDim;
        unsigned int ydim = YDim;
        mpREAL result[ydim], estabs[ydim];

        mpREAL xmin[xdim], xmax[xdim];
        for(int i=0; i<xdim; i++) {
            xmin[i] = 0;
            xmax[i] = 1;
        }
        LastState = 0;
        NEval = 0;
        nNAN = 0;
        
        size_t _MinPTS_, _RunPTS_;
        if(tn==0) {
            _RunPTS_ = RunPTS;
            MaxPTS = RunPTS * RunMAX;
            _MinPTS_ = MinPTS>0 ? MinPTS : _RunPTS_/10;
        } else {
            MaxPTS = tn;
            if(MaxPTS<10000) MaxPTS = 10000;
            _RunPTS_ = MaxPTS/5;
            _MinPTS_ = MinPTS>0 ? MinPTS : _RunPTS_/10;
        }
        StartTimer = time(NULL);
        StartTimer = time(NULL);

        Lib3_HCubatureMP::CPUCORES = omp_get_num_procs();
        int nok = Lib3_HCubatureMP::hcubature_v(ydim, Wrapper, this, xdim, xmin, xmax, _MinPTS_, _RunPTS_, MaxPTS, EpsAbs, EpsRel, result, estabs, tn==0 ? PrintHooker : NULL);

        if(nok) {
            mpREAL abs_res = 0;
            for(int j=0; j<ydim; j++) abs_res += result[j]*result[j];
            abs_res = sqrt(abs_res);
            mpREAL abs_est = 0;
            for(int j=0; j<ydim; j++) abs_est += estabs[j]*estabs[j];
            abs_est = sqrt(abs_est);
            mpREAL mpfr_eps = 10*mpfr::machine_epsilon();
            if( (abs_res < mpfr_eps) && (abs_est < mpfr_eps) ) {
                cout << ErrColor << "HCubatureMP Failed with 0 result returned!" << RESET << endl;
                return NaN;
            }
        }

        if(LastState==-1 && use_last) {
            for(int j=0; j<ydim; j++) {
                result[j] = LastResult[j];
                estabs[j] = LastAbsErr[j];
            }
            NEval = lastNRUN;
            nNAN = lastnNAN;
        }
        
        ex FResult = 0;
        bool any_NAN = false;
        for(int j=0; j<ydim; j++) {
            if(isnan(result[j])) {
                any_NAN = true;
                break;
            }
        }
        if(any_NAN) FResult = NaN;
        else {
            try{
                if(ydim==2) {
                    FResult = VE(mp2ex(result[0]), mp2ex(estabs[0])) + VE(mp2ex(result[1]), mp2ex(estabs[1])) * I;
                } else {
                    lst res;
                    for(int j=0; j<ydim; j++) res.append(VE(mp2ex(result[j]), mp2ex(estabs[j])));
                    FResult = res;
                }
            } catch(...) {
                FResult = NaN;
            }
        }
        
        // Check nNAN / NEval
        if(nNAN * 1000 > NEval && NEval>0) {
            cout << ErrColor << "NAN=" << nNAN << " v.s. RUN=" << NEval << RESET << endl;
        }
        
        return FResult;
    }

}
