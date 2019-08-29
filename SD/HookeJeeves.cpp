
#include "SD.h"
#include <math.h>

namespace HepLib {

dREAL HookeJeeves::best_nearby(dREAL* delta, dREAL* point, dREAL prevbest, int nvars) {
    dREAL z[nvars];
    dREAL minf, ftmp;
    int i;
    minf = prevbest;
    for (i = 0; i < nvars; i++) z[i] = point[i];
    for (i = 0; i < nvars; i++) {
        if(Exit) return (minf);
        z[i] = point[i] + delta[i];
        ftmp = ObjectWrapper(nvars, z);
        if (ftmp < minf) minf = ftmp;
        else {
            delta[i] = 0.0L - delta[i];
            z[i] = point[i] + delta[i];
            ftmp = ObjectWrapper(nvars, z);
            if (ftmp < minf) minf = ftmp;
            else z[i] = point[i];
        }
    }
    for (i = 0; i < nvars; i++) point[i] = z[i];
    return (minf);
}

int HookeJeeves::hooke(int nvars, dREAL* startpt, dREAL* endpt, dREAL rho, dREAL epsilon, int itermax) {
    dREAL delta[nvars];
    dREAL newf, fbefore, steplength, tmp;
    dREAL xbefore[nvars], newx[nvars];
    int i, j, keep;
    int iters, iadj;
    for (i = 0; i < nvars; i++) {
        newx[i] = xbefore[i] = startpt[i];
        delta[i] = std::fabs(startpt[i] * rho);
        if (delta[i] == 0.0) delta[i] = rho;
    }
    iadj = 0;
    steplength = rho;
    iters = 0;
    fbefore = ObjectWrapper(nvars, newx);
    newf = fbefore;
    while ((iters < itermax) && (steplength > epsilon)) {
        if(Exit) return (iters);
        iters++;
        iadj++;
        /* find best new point, one coord at a time */
        for (i = 0; i < nvars; i++) {
            newx[i] = xbefore[i];
        }
        newf = best_nearby(delta, newx, fbefore, nvars);
        /* if we made some improvements, pursue that direction */
        keep = 1;
        while ((newf < fbefore) && (keep == 1)) {
            if(Exit) return (iters);
            iadj = 0;
            for (i = 0; i < nvars; i++) {
                /* firstly, arrange the sign of delta[] */
                if (newx[i] <= xbefore[i]) delta[i] = 0.0 - std::fabs(delta[i]);
                else delta[i] = std::fabs(delta[i]);
                /* now, move further in this direction */
                tmp = xbefore[i];
                xbefore[i] = newx[i];
                newx[i] = newx[i] + newx[i] - tmp;
            }
            fbefore = newf;
            newf = best_nearby(delta, newx, fbefore, nvars);
            /* if the further (optimistic) move was bad.... */
            if (newf >= fbefore) break;
            /* make sure that the differences between the new */
            /* and the old points are due to actual */
            /* displacements; beware of roundoff errors that */
            /* might cause newf < fbefore */
            keep = 0;
            for (i = 0; i < nvars; i++) {
                keep = 1;
                if (std::fabs(newx[i] - xbefore[i]) > (0.5 * std::fabs(delta[i]))) break;
                else keep = 0;
            }
        }
        if ((steplength >= epsilon) && (newf >= fbefore)) {
            steplength = steplength * rho;
            for (i = 0; i < nvars; i++) {
                delta[i] *= rho;
            }
        }
    }
    for (i = 0; i < nvars; i++) endpt[i] = xbefore[i];
    return (iters);
}

dREAL HookeJeeves::ObjectWrapper(int nvars, dREAL* x) {
    for(int i=0; i<nvars; i++) {
        if(x[i]<LowerBound[i] || x[i]>UpperBound[i]) return 1E5;
    }
    return ObjectFunction(nvars, x, PL, LAS);
}

dREAL HookeJeeves::FindMinimum(int nvars, FunctionType func, dREAL *pl, dREAL *las, dREAL *UB, dREAL *LB, dREAL *IP, bool compare0) {
    ObjectFunction = func;
    PL = pl;
    LAS = las;
    
    if(UB != NULL) for(int i=0; i<nvars; i++) UpperBound[i] = UB[i];
    else for(int i=0; i<nvars; i++) UpperBound[i] = 1;
    
    if(LB != NULL) for(int i=0; i<nvars; i++) LowerBound[i] = LB[i];
    else for(int i=0; i<nvars; i++) LowerBound[i] = 0;
    
    dREAL RhoParameter = 0.95;
    dREAL EpsParameter = 1E-4;
    long long MaxParameter = 100000;
    
    #define savePTS 10
    #define tryPTS 2
    double mPoints[savePTS][nvars], mValue[savePTS];
    for(int i=0; i<savePTS; i++) mValue[i] = 1E5;
    int max_index = 0;
    dREAL pts[tryPTS+1];
    pts[0] = 1E-4;
    pts[tryPTS] = 1-1E-4;
    for(int i=1; i<tryPTS; i++) pts[i] = i*1.0/tryPTS;
    for(long long ii=0; ii<std::pow(tryPTS+1, nvars); ii++) {
        dREAL iPoints[nvars];
        int li = ii;
        for(int i=0; i<nvars; i++) {
            int mi = li % (1+tryPTS);
            iPoints[i] = LowerBound[i] + pts[mi]*(UpperBound[i]-LowerBound[i]);
            li /= (1+tryPTS);
        }
        auto tmp = ObjectWrapper(nvars, iPoints);
        if(mValue[max_index] > tmp) {
            mValue[max_index] = tmp;
            for(int j=0; j<nvars; j++) mPoints[max_index][j] = iPoints[j];
            max_index = 0;
            for(int j=0; j<savePTS && j<std::pow(tryPTS+1, nvars); j++) {
                if(mValue[j] > mValue[max_index]) max_index = j;
            }
        }
    }
    
    double ret = 1E5;
    for(int ii=0; ii<savePTS && ii<std::pow(tryPTS+1, nvars); ii++) {
        dREAL iPoints[nvars], oPoints[nvars];
        for(int i=0; i<nvars; i++) iPoints[i] = mPoints[ii][i];
        hooke(nvars, iPoints, oPoints, RhoParameter, EpsParameter, MaxParameter);
        auto tmp_ret = ObjectWrapper(nvars, oPoints);
        if(tmp_ret < ret) ret = tmp_ret;
    }
    
    return ret;
}

void HookeJeeves::ForceStop() {
    Exit = true;
}

void HookeJeeves::Minimize(int nvars, FunctionType func, dREAL *ip) {
    ObjectFunction = func;
    PL = NULL;
    LAS = NULL;
    
    for(int i=0; i<nvars; i++) UpperBound[i] = 100;
    for(int i=0; i<nvars; i++) LowerBound[i] = 0;
    
    dREAL EpsParameter = 1E-3;
    long long MaxParameter = 100000;
    
    dREAL iPoints[nvars], oPoints[nvars];
    for(int i=0; i<nvars; i++) iPoints[i] = ip[i];
    hooke(nvars, iPoints, oPoints, ILWrapper::hjRHO, EpsParameter, MaxParameter);
}


}
