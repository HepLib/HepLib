/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "AMF.h"
#include <cmath>

namespace HepLib {

    ex Fit::PolyFit(const lst & xs, const lst & ys, int lp) {
        unsigned int n = xs.nops();
        if(ys.nops() != n) throw Error("PolyFit: the size of xs is not the same as ys.");
        matrix X(n,n);
        for(int c=0; c<n; c++) {
            ex xp = 1;
            for(int r=0; r<n; r++) {
                X(r,c) = xp;
                xp *= xs.op(c);
            }
        }
        
        int nc = 1;
        bool is_y_lst = is_a<lst>(ys.op(0));
        if(is_y_lst) nc = ys.op(0).nops();
        matrix Y(n,nc);
        lst res;
        for(int c=0; c<nc; c++) {
            res.append(0);
            for(int r=0; r<n; r++) {
                Y(r,c) = is_y_lst ? ys.op(r).op(c) : ys.op(r);
                if(lp!=0) Y(r,c) /= pow(xs[r], lp);
            }
        }
        auto cmat = X.mul(X.transpose()).inverse(solve_algo::gauss).mul(X).mul(Y);
        
        for(int c=0; c<nc; c++) {
            for(int r=0; r<n; r++) res[c] += pow(ep,lp+r) * cmat(r,c);
        }
        return is_y_lst ? res : res.op(0);
    }
    
    // copy from Amflow.m
    ex Fit::GenerateNumericalConfig(int goal, int order) {
        // ep^number ~ ep^(-2 loop) 2^(-xorder)
        int loop = Internal.nops();
        int number = ceil(5*order/2.0+2*loop);
        if(number>100) {
            cout << "GenerateNumericalConfig: the given order is too large to evaluate." << endl;
            abort();
        }
        ex eps0 = pow(10, -1/ex(2)*loop-goal/ex(order+1));
        eps0 = Rationalize(NN(eps0, 3), 2);
        lst epslist;
        for(int i=0; i<number; i++) epslist.append(eps0*(1+(i+1)/ex(100)));
        int singlepre = ceil((number+2*loop)*(1/2.0*loop+goal/(order+1.0)));
        if(singlepre<30) singlepre = 30;
        int workpre = 2*singlepre;
        int xorder = 4*singlepre;
        return lst{epslist, workpre, xorder, eps0};
    }
    
    void Fit::operator()(int goal, int order, int rr) {
        auto config = GenerateNumericalConfig(goal, order);
        auto neps = config.op(0);
        int dp = ex2int(config.op(1));
        int xn = ex2int(config.op(2));
        auto nis = operator()(dp, xn, rr, ex_to<lst>(neps));
        int leading = -2*Internal.nops();
        lst eps = ex_to<lst>(neps);
        NIntegral = ex_to<lst>(PolyFit(eps, nis, leading));
        for(int i=0; i<NIntegral.nops(); i++) {
            NIntegral[i] = series_ex(NIntegral.op(i), ep, order+leading);
        }
    }
    
    lst Fit::operator()(int dp, int xn, int rr, lst eps) {
    
        if(!In_GiNaC_Parallel && Verbose>0) {
            cout << WHITE << "\\--Fit @ Start " << now() << RESET << endl;
            cout << "  \\--Configure: x^" << xn << ", working prec: " << dp << ", r=1/" << rr << " , ep numbers: " << eps.nops() << ", |ep| ~ " << eps.op(0) << endl;
            cout << endl;
        }
    
        lst nis;
        int cnt = 0, tot = eps.nops();
        for(auto nep : eps) {
            if(!In_GiNaC_Parallel && Verbose>0) {
                cnt++;
                cout << "  \\--[" << cnt << " / " << tot << "] @ eps = " << nep << endl;
            }
            ex nd = 4-2*nep;
            if(Mode==Modes::Single) {
                Single amf(nd);
            
                amf.Integral = Integral;
                amf.Propagator = Propagator;
                amf.Replacement = Replacement;
                amf.Internal = Internal;
                amf.External = External;
                
                amf.xn = xn;
                amf.dp = dp;
                amf.rr = rr;
                
                amf();
                nis.append(amf.NIntegral);
            } else if(Mode==Modes::Loop) {
                Loop amf(nd);
            
                amf.Integral = Integral;
                amf.Propagator = Propagator;
                amf.Replacement = Replacement;
                amf.Internal = Internal;
                amf.External = External;
                
                amf.xn = xn;
                amf.dp = dp;
                amf.rr = rr;
                
                amf();
                nis.append(amf.NIntegral);
            } else if(Mode==Modes::All) {
                All amf(nd);
            
                amf.Integral = Integral;
                amf.Propagator = Propagator;
                amf.Replacement = Replacement;
                amf.Internal = Internal;
                amf.External = External;
                
                amf.xn = xn;
                amf.dp = dp;
                amf.rr = rr;
                
                amf();
                nis.append(amf.NIntegral);
            } else throw Error("Not supported Mode.");
            
            if(!In_GiNaC_Parallel && Verbose>0) cout << endl;
        }
        
        if(!In_GiNaC_Parallel && Verbose>0) cout << WHITE << "\\--Fit @ Done " << now() << RESET << endl;
        return nis;
    }
    
    void Fit::Parallel(int goal, int order, int rr) {
        auto config = GenerateNumericalConfig(goal, order);
        auto neps = config.op(0);
        int dp = ex2int(config.op(1));
        int xn = ex2int(config.op(2));
        auto nis = Parallel(dp, xn, rr, ex_to<lst>(neps));
 
        int leading = -2*Internal.nops();
        lst eps = ex_to<lst>(neps);
        NIntegral = ex_to<lst>(PolyFit(eps, nis, leading));
        for(int i=0; i<NIntegral.nops(); i++) {
            NIntegral[i] = series_ex(NIntegral.op(i), ep, order+leading);
        }
    }
    
    lst Fit::Parallel(int dp, int xn, int rr, lst eps) {
    
        if(!In_GiNaC_Parallel && Verbose>0) {
            cout << WHITE << "\\--Fit @ Start " << now() << RESET << endl;
            cout << "  \\--Configure: x^" << xn << ", working prec: " << dp << ", r=1/" << rr << ", ep numbers: " << eps.nops() << ", |ep| ~ " << eps.op(0) << endl;
        }
    
        if(GiNaC_Parallel_Level>0) GiNaC_Parallel_Process = 0;
        auto res_vec = GiNaC_Parallel(eps.nops(), [&](int idx)->ex {
            auto nep = eps.op(idx);
            ex nd = 4-2*nep;
            if(Mode==Modes::Single) {
                Single amf(nd);
            
                amf.Integral = Integral;
                amf.Propagator = Propagator;
                amf.Replacement = Replacement;
                amf.Internal = Internal;
                amf.External = External;
                
                amf.xn = xn;
                amf.dp = dp;
                amf.rr = rr;
                
                amf();
                return amf.NIntegral;
            } else if(Mode==Modes::Loop) {
                Loop amf(nd);
            
                amf.Integral = Integral;
                amf.Propagator = Propagator;
                amf.Replacement = Replacement;
                amf.Internal = Internal;
                amf.External = External;
                
                amf.xn = xn;
                amf.dp = dp;
                amf.rr = rr;
                
                amf();
                return amf.NIntegral;
            } else if(Mode==Modes::All) {
                All amf(nd);
            
                amf.Integral = Integral;
                amf.Propagator = Propagator;
                amf.Replacement = Replacement;
                amf.Internal = Internal;
                amf.External = External;
                
                amf.xn = xn;
                amf.dp = dp;
                amf.rr = rr;
                
                amf();
                return amf.NIntegral;
            } else throw Error("Not supported Mode.");
        }, "Fit");
        
        if(!In_GiNaC_Parallel && Verbose>0) cout << WHITE << "\\--Fit @ Done " << now() << RESET << endl;
        
        lst nis;
        for(auto item : res_vec) nis.append(item);

        return nis;
    }
    
}

