/**
 * @file
 * @brief Functions for Renormalization Constants
 */
 
#include "QCD.h"

namespace HepLib::QCD {

    //-----------------------------------------------------------
    // FF namespace for Fragmentation Function
    //-----------------------------------------------------------
    namespace FF {
        
        /**
         * @brief the prefactor for the z-Integration
         * arXiv:1208.5301v2 [hep-ph] 9 Nov 2012
         * PHYSICAL REVIEW D 91, 074013 (2015)
         * @param mode 1-Gluon, 2-Quark
         * @param tls the number of particles in cut
         * @param SF the symmetry factor n, will used as 1/n
         * @return the prefactor for the z-Integration
         */
        ex zIntFactor(int mode, int tls, const ex SF) {
            Symbol kp("kp"), z("z"), M("M");
            ex ret = 0, NCS = 1;
            if(mode==1) NCS = -pow(z,1-2*ep)/(16*(2-2*ep)*kp*Pi);
            if(mode==2) NCS = pow(z,(1-2*ep))/(8*3*Pi);
            else throw Error("zIntFactor: Not supported mode "+to_string(mode));
            if(tls>0) ret = NCS*M/SF*(-((pow(2,(2-tls-(3-2*ep)*tls))*pow(Pi,(1-(3-2*ep)*tls)))/(kp*(-1+z))));
            else ret = NCS*M*4*Pi/kp;
            return Symbol::set_all(ret);
        }
        
        /**
         * @brief dPhi (c1 p.k + c0)^n, NOT including the zIntFactor yet
         * @param c1 coefficient of p.k
         * @param c0 constant term
         * @param n the exponent
         * @return the z integration
         */
        ex zIntegrate(const ex & c1, const ex & c0, const ex & n, const ex k2, const ex& p2) {
            Symbol z("z"), kp("kp"), pp("pp");
            ex ret =  pow(2,n)*pow(Pi,1-ep)*pow(2*c0+c1*kp*(p2-p2*z)*pow(pp,-1)+c1*k2*pp*pow(kp-kp*z,-1),-n)*pow(c1*pow(pp,2)*pow(-2*c0*kp*pp*(-1+z)+c1*(k2*pow(pp,2)+p2*pow(kp,2)*pow(-1+z,2)),-1),-1+ep)*pow(tgamma(n),-1)*tgamma(-1+ep+n);
            return Symbol::set_all(ret);
        }
        
        ex _Anti5(const ex & expr) { // for internal usage only
            if(!DGamma::has(expr)) return expr;
            if(is_a<DGamma>(expr)) {
                DGamma g = ex_to<DGamma>(expr);
                if(is_a<Vector>(g.pi) || is_a<Index>(g.pi)) return -expr;
                else if(is_zero(g.pi-1) || is_zero(g.pi-5)) return expr;
            }
            
            if(is_a<add>(expr)) {
                ex ret = 0;
                for(auto item : expr) ret += _Anti5(item);
                return ret;
            } else if(is_a<mul>(expr) || is_a<ncmul>(expr)) {
                ex ret = 1;
                for(auto item : expr) ret *= _Anti5(item);
                return ret;
            } 
            cout << DGamma::has(expr) << " : " << expr << endl;
            throw Error("_Anti5: unexpected region.");
            return 0;
        }
        
        ex Anti5R(const ex & expr) {
            static ex g5 = GAS(5);
            static MapFunction anti5([](const ex & e, MapFunction self)->ex{
                if(!e.has(g5)) return e;
                else if(is_a<ncmul>(e)) {
                    ex eL=1, eR=1;
                    bool found = false;
                    for(auto item : e) {
                        if(!found && item.is_equal(g5)) {
                            found = true;
                            continue;
                        }
                        if(found) eR = eR * item;
                        else eL = eL * item;
                    }
                    eR = _Anti5(eR);
                    if(eR.has(GAS(5))) eR = Anti5R(eR);
                    return eL * eR * g5;
                } else return e.map(self);
                return 0;
            });
            return anti5(expr);
        }
        
    }
}

