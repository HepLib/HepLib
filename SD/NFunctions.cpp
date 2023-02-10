/**
 * @file
 * @brief Functions used in generated C++ numerical code
 */

#include "NFunctions.h"

int NRCLog = 5;

const dREAL dPi = 3.1415926535897932384626433832795028841971693993751L;
const dREAL dEuler = 0.57721566490153286060651209008240243104215933593992L;
const dCOMPLEX diEpsilon = complex<dREAL>(0, 100*LDBL_EPSILON);

const qREAL qPi = 3.1415926535897932384626433832795028841971693993751Q;
const qREAL qEuler = 0.57721566490153286060651209008240243104215933593992Q;
const qCOMPLEX qiEpsilon = 100.Qi*FLT128_EPSILON;

mpREAL mpPi;
mpREAL mpEuler;
mpCOMPLEX mpiEpsilon;

void X2Z(int nfxs, dREAL(*f)(const dREAL*,const dREAL*), dREAL(*Df)(const int,const dREAL*,const dREAL*),
    const dREAL* x, dCOMPLEX* z, dCOMPLEX* r, dREAL* dff, const dREAL* pl, const dREAL* las) {
    dCOMPLEX ilas[nfxs];
    for(int i=0; i<nfxs; i++) ilas[i] = complex<dREAL>(0.L, las[i]);
    dff[nfxs] = f(x,pl);
    for(int i=0; i<nfxs; i++) dff[i] = Df(i,x,pl);
    for(int i=0; i<nfxs; i++) r[i] = dff[i]*ilas[i];
    for(int i=0; i<nfxs; i++) z[i] = x[i]-x[i]*(1.L-x[i])*r[i];
}

void Mat(int nfxs, dREAL(*DDf)(const int,const int,const dREAL*,const dREAL*), 
    dCOMPLEX* mat, const dREAL* x, const dREAL* dff, const dREAL* pl, const dREAL* las) {
    dCOMPLEX ilas[nfxs];
    for(int i=0; i<nfxs; i++) ilas[i] = complex<long double>(0.L, las[i]);
    dREAL ddf[nfxs][nfxs];
    for(int i=0; i<nfxs; i++) {
        for(int j=0; j<nfxs; j++) ddf[i][j] = DDf(i,j,x,pl);
    }
    for(int i=0; i<nfxs; i++) {
        for(int j=0; j<nfxs; j++) {
            int ij = i*nfxs+j;
            if(i!=j) mat[ij] = 0;
            else mat[ij] = 1.L-(1.L-2.L*x[i])*dff[i]*ilas[i];
            mat[ij] = mat[ij]-x[i]*(1.L-x[i])*ddf[i][j]*ilas[i];
        }
    }
}

dCOMPLEX MatDet(dCOMPLEX mat[], int n) {
    dREAL epsilon = LDBL_EPSILON * 100;
    bool is_zero = false;
    int s=1;
    for(int i=0; i<n-1; i++) {
        if(fabs(mat[i*n+i])<epsilon) {
            bool is_zero = true;
            for(int j=i+1; j<n; j++) {
                if(fabs(mat[i*n+j])>epsilon) {
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
    dCOMPLEX ret = log(xs[n]);
    if(n<2) return ret;
    auto eps = 100*LDBL_EPSILON;
    int total=0;
    int ReIm[n+1][2];
    for(int k=0; k<=n; k++) {
        auto curR = xs[k].real();
        auto curI = xs[k].imag();
        auto absR = fabs(curR);
        auto absI = fabs(curI);
        if(absR<10*eps*absI || absI<10*eps*absR) continue; 
        ReIm[total][0] = curR>0 ? 1 : -1;
        ReIm[total][1] = curI>0 ? 1 : -1;
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

void X2Z(int nfxs, qREAL(*f)(const qREAL*,const qREAL*), qREAL(*Df)(const int,const qREAL*,const qREAL*),
    const qREAL* x, qCOMPLEX* z, qCOMPLEX* r, qREAL* dff, const qREAL* pl, const qREAL* las) {
    qCOMPLEX ilas[nfxs];
    for(int i=0; i<nfxs; i++) ilas[i] = las[i] * 1.Qi;
    dff[nfxs] = f(x,pl);
    for(int i=0; i<nfxs; i++) dff[i] = Df(i,x,pl);
    for(int i=0; i<nfxs; i++) r[i] = dff[i]*ilas[i];
    for(int i=0; i<nfxs; i++) z[i] = x[i]-x[i]*(1.Q-x[i])*r[i];
}

void Mat(int nfxs, qREAL(*DDf)(const int,const int,const qREAL*,const qREAL*), 
    qCOMPLEX *mat, const qREAL* x, const qREAL* dff, const qREAL *pl, const qREAL *las) {
    qCOMPLEX ilas[nfxs];
    for(int i=0; i<nfxs; i++) ilas[i] = las[i] * 1.Qi;
    qREAL ddf[nfxs][nfxs];
    for(int i=0; i<nfxs; i++) {
        for(int j=0; j<nfxs; j++) ddf[i][j] = DDf(i,j,x,pl);
    }
    for(int i=0; i<nfxs; i++) {
        for(int j=0; j<nfxs; j++) {
            int ij = i*nfxs+j;
            if(i!=j) mat[ij] = 0;
            else mat[ij] = 1.Q-(1.Q-2.Q*x[i])*dff[i]*ilas[i];
            mat[ij] = mat[ij]-x[i]*(1.Q-x[i])*ddf[i][j]*ilas[i];
        }
    }
}

qCOMPLEX MatDet(qCOMPLEX mat[], int n) {
    qREAL epsilon = FLT128_EPSILON * 100;
    bool is_zero = false;
    int s=1;
    for(int i=0; i<n-1; i++) {
        if(cabsq(mat[i*n+i])<epsilon) {
            bool is_zero = true;
            for(int j=i+1; j<n; j++) {
                if(cabsq(mat[i*n+j])>epsilon) {
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
    qCOMPLEX ret = log(xs[n]);
    if(n<2) return ret;
    auto eps = 100*FLT128_EPSILON;
    int total=0;
    int ReIm[n+1][2];
    for(int k=0; k<=n; k++) {
        auto curR = crealq(xs[k]);
        auto curI = cimagq(xs[k]);
        auto absR = fabsq(curR);
        auto absI = fabsq(curI);
        if(absR<10*eps*absI || absI<10*eps*absR) continue; 
        ReIm[total][0] = curR>0 ? 1 : -1;
        ReIm[total][1] = curI>0 ? 1 : -1;
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

void X2Z(int nfxs, mpREAL(*f)(const mpREAL*,const mpREAL*), mpREAL(*Df)(const int,const mpREAL*,const mpREAL*),
    const mpREAL* x, mpCOMPLEX* z, mpCOMPLEX* r, mpREAL* dff, const mpREAL* pl, const mpREAL* las) {
    mpCOMPLEX ilas[nfxs];
    for(int i=0; i<nfxs; i++) ilas[i] = complex<mpREAL>(mpREAL(0), las[i]);
    dff[nfxs] = f(x,pl);
    for(int i=0; i<nfxs; i++) dff[i] = Df(i,x,pl);
    for(int i=0; i<nfxs; i++) r[i] = dff[i]*ilas[i];
    for(int i=0; i<nfxs; i++) z[i] = x[i]-x[i]*(1-x[i])*r[i];
}

void Mat(int nfxs, mpREAL(*DDf)(const int, const int,const mpREAL*,const mpREAL*), 
    mpCOMPLEX *mat, const mpREAL* x, const mpREAL* dff, const mpREAL *pl, const mpREAL *las) {
    mpCOMPLEX ilas[nfxs];
    for(int i=0; i<nfxs; i++) ilas[i] = complex<mpREAL>(mpREAL(0), las[i]);
    mpREAL ddf[nfxs][nfxs];
    for(int i=0; i<nfxs; i++) {
        for(int j=0; j<nfxs; j++) ddf[i][j] = DDf(i,j,x,pl);
    }
    for(int i=0; i<nfxs; i++) {
        for(int j=0; j<nfxs; j++) {
            int ij = i*nfxs+j;
            if(i!=j) mat[ij] = 0;
            else mat[ij] = mpREAL(1)-(1-2*x[i])*dff[i]*ilas[i];
            mat[ij] = mat[ij]-x[i]*(1-x[i])*ddf[i][j]*ilas[i];
        }
    }
}

mpCOMPLEX MatDet(mpCOMPLEX mat[], int n) {
    bool is_zero = false;
    mpREAL epsilon = mpfr::machine_epsilon() * 100;
    int s=1;
    for(int i=0; i<n-1; i++) {
        if(abs(mat[i*n+i])<epsilon) {
            bool is_zero = true;
            for(int j=i+1; j<n; j++) {
                if(abs(mat[i*n+j])>epsilon) {
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
    mpCOMPLEX ret = log(xs[n]);
    if(n<2) return ret;
    auto eps = mpfr::machine_epsilon() * 100;
    int total=0;
    int ReIm[n+1][2];
    for(int k=0; k<=n; k++) {
        auto curR = xs[k].real();
        auto curI = xs[k].imag();
        auto absR = abs(curR);
        auto absI = abs(curI);
        if(absR<10*eps*absI || absI<10*eps*absR) continue; 
        ReIm[total][0] = curR>0 ? 1 : -1;
        ReIm[total][1] = curI>0 ? 1 : -1;
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

