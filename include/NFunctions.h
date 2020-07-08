 
#pragma once

#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <iostream>
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

extern int NRCLog;
extern const dREAL dPi;
extern const dREAL dEuler;
extern const dCOMPLEX diEpsilon;
extern const qREAL qPi;
extern const qREAL qEuler;
extern const qCOMPLEX qiEpsilon;
extern mpREAL mpPi;
extern mpREAL mpEuler;
extern mpCOMPLEX mpiEpsilon;

void X2Z(int nfxs, dREAL(*f)(const dREAL*,const dREAL*), dREAL(*Df)(const int,const dREAL*,const dREAL*),
    const dREAL* x, dCOMPLEX* z, dCOMPLEX* r, dREAL* dff, const dREAL* pl, const dREAL* las);
void Mat(int nfxs, dREAL(*DDf)(const int,const int,const dREAL*,const dREAL*), 
    dCOMPLEX* mat, const dREAL* x, const dREAL* dff, const dREAL* pl, const dREAL* las);
    
void X2Z(int nfxs, qREAL(*f)(const qREAL*,const qREAL*), qREAL(*Df)(const int,const qREAL*,const qREAL*),
    const qREAL* x, qCOMPLEX* z, qCOMPLEX* r, qREAL* dff, const qREAL* pl, const qREAL* las);
void Mat(int nfxs, qREAL(*DDf)(const int,const int,const qREAL*,const qREAL*), 
    qCOMPLEX *mat, const qREAL* x, const qREAL* dff, const qREAL *pl, const qREAL *las);
    
void X2Z(int nfxs, mpREAL(*f)(const mpREAL*,const mpREAL*), mpREAL(*Df)(const int,const mpREAL*,const mpREAL*),
    const mpREAL* x, mpCOMPLEX* z, mpCOMPLEX* r, mpREAL* dff, const mpREAL* pl, const mpREAL* las);
void Mat(int nfxs, mpREAL(*DDf)(const int, const int,const mpREAL*,const mpREAL*), 
    mpCOMPLEX *mat, const mpREAL* x, const mpREAL* dff, const mpREAL *pl, const mpREAL *las);
    
dCOMPLEX MatDet(dCOMPLEX mat[], int n);
qCOMPLEX MatDet(qCOMPLEX mat[], int n);
mpCOMPLEX MatDet(mpCOMPLEX mat[], int n);

dCOMPLEX RCLog(dCOMPLEX xs[], int n);
qCOMPLEX RCLog(qCOMPLEX xs[], int n);
mpCOMPLEX RCLog(mpCOMPLEX xs[], int n);

inline dREAL expt(dREAL a, dREAL b) { return pow(a,b); }
inline dCOMPLEX expt(dCOMPLEX a, dREAL b) { return pow(a,b); }
inline dREAL recip(dREAL a) { return 1.L/a; }
inline dCOMPLEX recip(dCOMPLEX a) { return 1.L/a; }

inline qREAL expt(qREAL a, qREAL b) { return powq(a,b); }
inline qCOMPLEX expt(qCOMPLEX a, qREAL b) { return cpowq(a,b); }
inline qREAL recip(qREAL a) { return 1.Q/a; }
inline qCOMPLEX recip(qCOMPLEX a) { return 1.Q/a; }

inline mpREAL expt(mpREAL a, mpREAL b) { return pow(a,b); }
inline mpCOMPLEX expt(mpCOMPLEX a, mpREAL b) { return pow(a,b); }
inline mpREAL recip(mpREAL a) { return mpREAL(1)/a; }
inline mpCOMPLEX recip(mpCOMPLEX a) { return mpREAL(1)/a; }

inline qREAL pow(qREAL x, qREAL y) { return powq(x, y); }
inline qREAL log(qREAL x) { return logq(x); }
inline qREAL exp(qREAL x) { return expq(x); }
inline qCOMPLEX pow(qCOMPLEX x, qREAL y) { return cpowq(x, y); }
inline qCOMPLEX log(qCOMPLEX x) { return clogq(x); }
inline qCOMPLEX exp(qCOMPLEX x) { return cexpq(x); }
