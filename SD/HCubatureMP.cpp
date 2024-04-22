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
            mpfr::mpreal::set_default_prec(mpfr::digits2bits(self->MPDigits));
            self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->mpParameter, self->mpLambda);
            // Final Check NaN/Inf
            bool ok = true;
            for(int j=0; j<ydim; j++) {
                mpREAL ytmp = y[i*ydim+j];
                if(isnan(ytmp) || isinf(ytmp)) { ok = false; break; }
            }
            if(!ok && (self->IntegrandMP!=NULL)) {
                mpfr_free_cache();
                mpfr::mpreal::set_default_prec(mpfr::digits2bits(self->MPDigits*10));
                self->IntegrandMP(xdim, x+i*xdim, ydim, y+i*ydim, self->mpParameter, self->mpLambda);
                mpfr::mpreal::set_default_prec(mpfr::digits2bits(self->MPDigits));
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
            if(self->ReIm == 1) y[i*ydim+1] = 0;
            else if(self->ReIm == 2) y[i*ydim+0] = 0;
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
            auto r0 = result[0];
            auto r1 = result[1];
            auto e0 = epsabs[0].toString(3);
            auto e1 = epsabs[1].toString(3);
            cout << "     N: " << (*nrun) << ", ";
            if(self->ReIm==3 || self->ReIm==1) cout << "[" << r0 << ", " << e0 << "]";
            if(self->ReIm==3 || self->ReIm==2) cout << "+I*[" << r1 << ", " << e1 << "]";
            cout << endl;
        }
        if((*nrun-self->NEval) >= self->RunPTS) self->NEval = *nrun;
        
        if((isnan(result[0]) || isnan(result[1]) || isnan(epsabs[0]) || isnan(epsabs[1])) || (isinf(result[0]) || isinf(result[1]) || isinf(epsabs[0]) || isinf(epsabs[1]))) {
            self->NEval = *nrun;
            *nrun = self->MaxPTS + 1979;
            if(self->LastState>0) self->LastState = -1;
            if(Verbose>10) cout << ErrColor << "     Exit with NaN, LastN=" << self->lastNRUN << RESET << endl;
            return;
        }
        
        if(epsabs[0] > 1E30*self->EpsAbs || epsabs[1] > 1E30*self->EpsAbs) {
            self->NEval = *nrun;
            *nrun = self->MaxPTS + 1979;
            if(self->LastState>0) self->LastState = -1;
            if(Verbose>10) cout << WarnColor << "     Exit with EpsAbs, LastN=" << self->lastNRUN << RESET << endl;
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
        if(rExit && iExit && (*nrun)>self->MinPTS) {
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
        if(mpfr_buildopt_tls_p()<=0) throw Error("Integrate: mpfr_buildopt_tls_p()<=0.");
        mpfr_free_cache();
        mpfr::mpreal::set_default_prec(mpfr::digits2bits(MPDigits));
        mpPi = mpfr::const_pi();
        mpEuler = mpfr::const_euler();
        mpiEpsilon = complex<mpREAL>(0,mpfr::machine_epsilon()*100);
        
        unsigned int xdim = XDim;
        unsigned int ydim = 2;
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
            mpREAL abs_res = sqrt(result[0]*result[0]+result[1]*result[1]);
            mpREAL abs_est = sqrt(estabs[0]*estabs[0]+estabs[1]*estabs[1]);
            mpREAL mpfr_eps = 10*mpfr::machine_epsilon();
            if( (abs_res < mpfr_eps) && (abs_est < mpfr_eps) ) {
                cout << ErrColor << "HCubatureMP Failed with 0 result returned!" << RESET << endl;
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
        
        // Check nNAN / NEval
        if(nNAN * 1000 > NEval && NEval>0) {
            cout << ErrColor << "NAN=" << nNAN << " v.s. RUN=" << NEval << RESET << endl;
        }
        
        return FResult;
    }

}
