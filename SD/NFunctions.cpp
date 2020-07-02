/**
 * @file
 * @brief Matrix Determinat, used in Numerical Integration
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
  
#include <math.h>
#include <complex>
extern "C" {
#include <quadmath.h>
}
#include "mpreal.h"

int NRCLog = 7;

using namespace std;
typedef __float128 qREAL;
typedef __complex128 qCOMPLEX;
typedef long double dREAL;
typedef complex<long double> dCOMPLEX;
typedef mpfr::mpreal mpREAL;
typedef complex<mpREAL> mpCOMPLEX;

dREAL expt(dREAL a, dREAL b) { return pow(a,b); }
dCOMPLEX expt(dCOMPLEX a, dREAL b) { return pow(a,b); }
dREAL recip(dREAL a) { return 1.L/a; }
dCOMPLEX recip(dCOMPLEX a) { return 1.L/a; }

qREAL expt(qREAL a, qREAL b) { return powq(a,b); }
qCOMPLEX expt(qCOMPLEX a, qREAL b) { return cpowq(a,b); }
qREAL recip(qREAL a) { return 1.Q/a; }
qCOMPLEX recip(qCOMPLEX a) { return 1.Q/a; }

mpREAL expt(mpREAL a, mpREAL b) { return pow(a,b); }
mpCOMPLEX expt(mpCOMPLEX a, mpREAL b) { return pow(a,b); }
mpREAL recip(mpREAL a) { return mpREAL(1)/a; }
mpCOMPLEX recip(mpCOMPLEX a) { return mpREAL(1)/a; }

qREAL pow(qREAL x, qREAL y) { return powq(x, y); }
qREAL log(qREAL x) { return logq(x); }
qCOMPLEX pow(qCOMPLEX x, qREAL y) { return cpowq(x, y); }
qCOMPLEX log(qCOMPLEX x) { return clogq(x); }
qCOMPLEX exp(qCOMPLEX x) { return cexpq(x); }

dREAL dPi = 3.1415926535897932384626433832795028841971693993751L;
dREAL dEuler = 0.57721566490153286060651209008240243104215933593992L;
dCOMPLEX diEpsilon = complex<dREAL>(0, 1.E-50);

dCOMPLEX MatDet(dCOMPLEX mat[], int n) {
    bool is_zero = false;
    int s=1;
    for(int i=0; i<n-1; i++) {
        if(fabs(mat[i*n+i])<1.0E-15) {
            bool is_zero = true;
            for(int j=i+1; j<n; j++) {
                if(fabs(mat[i*n+j])>1.0E-15) {
                    for(int k=0; k<n; k++) {
                        auto tmp = mat[k*n+j];
                        mat[k*n+j] = mat[k*n+i];
                        mat[k*n+i] = tmp;
                    }
                    is_zero = false;
                    s=-s;
                    break;
                }
            }
            if(is_zero) return 0;
        }
        for(int k=i+1; k<n; k++) {
            auto m = mat[k*n+i]/mat[i*n+i];
            for(int j=0; j<n; j++) mat[k*n+j] = mat[k*n+j] - m*mat[i*n+j];
        }
    }
    dCOMPLEX ret = s;
    for(int k=0; k<n; k++) ret *= mat[k*n+k];
    return ret;
}

dCOMPLEX RCLog(dCOMPLEX xs[], int n) {
    auto eps = LDBL_EPSILON;
    dCOMPLEX ret = log(xs[n-1]);
    if(n<2) return ret;
    int total=0;
    int ReIm[n][2];
    for(int k=0; k<n; k++) {
        auto curR = xs[k].real();
        auto curI = xs[k].imag();
        auto absR = fabs(curR);
        auto absI = fabs(curI);
        if(absI<10*eps*absR && k==0 && absR<0) ReIm[total][1] = -1;
        else if(absR<10*eps*absI || absI<10*eps*absR) continue; 
        else ReIm[total][1] = curI>0 ? 1 : -1;
        ReIm[total][0] = curR>0 ? 1 : -1;
        total++;
    }
    
    int cutN = 0;
    for(int k=0; k<total-1; k++) {
        if(ReIm[k][0]*ReIm[k+1][0]<0 && ReIm[k][1]*ReIm[k+1][1]<0) return nanl("");
        if(ReIm[k][0]<0 && ReIm[k+1][0]<0 && ReIm[k][1]*ReIm[k+1][1]<0) {
            if(ReIm[k][1]>0) cutN++;
            else cutN--;
        }
    }
    if(cutN!=0) ret += complex<dREAL>(0,cutN * 2 * dPi);
    return ret;
}

qREAL qPi = 3.1415926535897932384626433832795028841971693993751Q;
qREAL qEuler = 0.57721566490153286060651209008240243104215933593992Q;
qCOMPLEX qiEpsilon = 1.E-50Qi;

qCOMPLEX MatDet(qCOMPLEX mat[], int n) {
    bool is_zero = false;
    int s=1;
    for(int i=0; i<n-1; i++) {
        if(cabsq(mat[i*n+i])<1.0E-25) {
            bool is_zero = true;
            for(int j=i+1; j<n; j++) {
                if(cabsq(mat[i*n+j])>1.0E-25) {
                    for(int k=0; k<n; k++) {
                        auto tmp = mat[k*n+j];
                        mat[k*n+j] = mat[k*n+i];
                        mat[k*n+i] = tmp;
                    }
                    is_zero = false;
                    s=-s;
                    break;
                }
            }
            if(is_zero) return 0;
        }
        for(int k=i+1; k<n; k++) {
            auto m = mat[k*n+i]/mat[i*n+i];
            for(int j=0; j<n; j++) mat[k*n+j] = mat[k*n+j] - m*mat[i*n+j];
        }
    }
    qCOMPLEX ret = s;
    for(int k=0; k<n; k++) ret *= mat[k*n+k];
    return ret;
}

qCOMPLEX RCLog(qCOMPLEX xs[], int n) {
    auto eps = FLT128_EPSILON;
    qCOMPLEX ret = log(xs[n-1]);
    if(n<2) return ret;
    int total=0;
    int ReIm[n][2];
    for(int k=0; k<n; k++) {
        auto curR = crealq(xs[k]);
        auto curI = cimagq(xs[k]);
        auto absR = fabsq(curR);
        auto absI = fabsq(curI);
        if(absI<10*eps*absR && k==0 && absR<0) ReIm[0][1] = -1;
        else if(absR<10*eps*absI || absI<10*eps*absR) continue; 
        else ReIm[total][1] = curI>0 ? 1 : -1;
        ReIm[total][0] = curR>0 ? 1 : -1;
        total++;
    }
    
    int cutN = 0;
    for(int k=0; k<total-1; k++) {
        if(ReIm[k][0]*ReIm[k+1][0]<0 && ReIm[k][1]*ReIm[k+1][1]<0) return nanl("");
        if(ReIm[k][0]<0 && ReIm[k+1][0]<0 && ReIm[k][1]*ReIm[k+1][1]<0) {
            if(ReIm[k][1]>0) cutN++;
            else cutN--;
        }
    }
    if(cutN!=0) ret += cutN * qPi * 2.Qi;
    return ret;
}

mpREAL mpPi;
mpREAL mpEuler;
mpCOMPLEX mpiEpsilon = complex<mpREAL>(0, 1.E-50);

mpCOMPLEX MatDet(mpCOMPLEX mat[], int n) {
    bool is_zero = false;
    int s=1;
    for(int i=0; i<n-1; i++) {
        if(abs(mat[i*n+i])<1.0E-35) {
            bool is_zero = true;
            for(int j=i+1; j<n; j++) {
                if(abs(mat[i*n+j])>1.0E-35) {
                    for(int k=0; k<n; k++) {
                        auto tmp = mat[k*n+j];
                        mat[k*n+j] = mat[k*n+i];
                        mat[k*n+i] = tmp;
                    }
                    is_zero = false;
                    s=-s;
                    break;
                }
            }
            if(is_zero) return mpREAL(0);
        }
        for(int k=i+1; k<n; k++) {
            auto m = mat[k*n+i]/mat[i*n+i];
            for(int j=0; j<n; j++) mat[k*n+j] = mat[k*n+j] - m*mat[i*n+j];
        }
    }
    mpCOMPLEX ret = mpREAL(s);
    for(int k=0; k<n; k++) ret *= mat[k*n+k];
    return ret;
}

mpCOMPLEX RCLog(mpCOMPLEX xs[], int n) {
    auto eps = mpfr::machine_epsilon(xs[0].real());
    mpCOMPLEX ret = log(xs[n-1]);
    if(n<2) return ret;
    int total=0;
    int ReIm[n][2];
    for(int k=0; k<n; k++) {
        auto curR = xs[k].real();
        auto curI = xs[k].imag();
        auto absR = abs(curR);
        auto absI = abs(curI);
        if(absI<10*eps*absR && k==0 && absR<0) ReIm[total][1] = -1;
        else if(absR<10*eps*absI || absI<10*eps*absR) continue; 
        else ReIm[total][1] = curI>0 ? 1 : -1;
        ReIm[total][0] = curR>0 ? 1 : -1;
        total++;
    }
    
    int cutN = 0;
    for(int k=0; k<total-1; k++) {
        if(ReIm[k][0]*ReIm[k+1][0]<0 && ReIm[k][1]*ReIm[k+1][1]<0) return mpREAL(nanl(""));
        if(ReIm[k][0]<0 && ReIm[k+1][0]<0 && ReIm[k][1]*ReIm[k+1][1]<0) {
            if(ReIm[k][1]>0) cutN++;
            else cutN--;
        }
    }
    if(cutN!=0) ret += complex<mpREAL>(0,cutN * 2 * mpPi);
    return ret;
}

