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

#define Pi 3.1415926535897932384626433832795028841971693993751L
#define Euler 0.57721566490153286060651209008240243104215933593992L

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

dCOMPLEX MatDetL(dCOMPLEX mat[], int n) {
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

#undef Pi
#undef Euler
#define Pi 3.1415926535897932384626433832795028841971693993751Q
#define Euler 0.57721566490153286060651209008240243104215933593992Q

qCOMPLEX MatDetQ(qCOMPLEX mat[], int n) {
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

#undef Pi
#undef Euler
#define Pi mpREAL("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068")
#define Euler mpREAL("0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495")

mpCOMPLEX MatDetMP(mpCOMPLEX mat[], int n) {
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
