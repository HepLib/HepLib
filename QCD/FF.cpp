/**
 * @file
 * @brief Functions for Renormalization Constants
 */
 
#include "QCD.h"
#include "SD.h"

namespace HepLib::QCD {

    //-----------------------------------------------------------
    // FF namespace for Fragmentation Function
    //-----------------------------------------------------------
    namespace FF {
        
        /**
         * @brief Feyman Parameterization, without zIntFactor
         * @param ls loop lst
         * @param tls tloop lst
         * @param ps propagator lst
         * @param ns exponent lst
         * @param lr loop replacement
         * @param tlr tLoop replacement
         * @param nr numeric replacement
         * @return Feyman Parameterization
         */
        ex FeynParameterize(const lst & ls, const lst & tls, const lst & ps, const lst & ns, const lst & lr, const lst & tlr, const lst & nr) {
        
            SD::FeynmanParameter fp;
            fp.LoopMomenta = ex_to<lst>(ls);
            fp.tLoopMomenta = ex_to<lst>(tls);
            fp.Propagator = ex_to<lst>(ps);
            fp.Exponent = ex_to<lst>(ns);
            
            for(auto vv : lr) fp.lReplacement[vv.op(0)] = vv.op(1);
            for(auto vv : tlr) fp.tReplacement[vv.op(0)] = vv.op(1);
            for(auto vv : nr) fp.nReplacement[vv.op(0)] = vv.op(1);
                        
            SD::SecDec work;
            work.Initialize(fp);
            
            ex ret = 0;
            for(auto fe : work.FunExp) {
                ret += WF(fe.op(0),fe.op(1),fe.op(3));
            }
            return ret;
        }
    
        /**
         * @brief Eikonal Propagator, left side w.r.t cutting line
         * @param e expression with head of Propagator
         * @param n the direction vector
         * @param mode 0 for gluon, others for quark/anti-quark
         * @return Eikonal Propagator
         */
        ex eikonalPropagator(ex e, ex n, int mode) {
            auto fi1 = e.op(0).op(1);
            auto fi2 = e.op(1).op(1);
            auto mom = e.op(2);
            if(mode==0) return I * SP(CI(fi1),CI(fi2)) / (SP(n,mom)+iEpsilon);
            else return I * GMat(GAS(1), DI(fi1),DI(fi2)) * SP(TI(fi1),TI(fi2)) / (SP(n,mom)+iEpsilon);
        }
        
        /**
         * @brief Eikonal Propagator, right side version w.r.t to cutting line
         * @param e expression with head of Propagator
         * @param n the direction vector
         * @param mode 0 for gluon, others for quark/anti-quark
         * @return Eikonal Propagator
         */
        ex eikonalPropagatorR(ex e, ex n, int mode) {
            auto fi1 = e.op(0).op(1);
            auto fi2 = e.op(1).op(1);
            auto mom = e.op(2);
            if(mode==0) return -I * SP(CI(fi1),CI(fi2)) / (SP(n,mom)-iEpsilon);
            else return -I * GMat(GAS(1), DI(fi1),DI(fi2)) * SP(TI(fi1),TI(fi2)) / (SP(n,mom)-iEpsilon);
        }
        
        /**
         * @brief Eikonal-Gloun Vertex, left side w.r.t. cutting line
         * @param e expression with head of Vertex
         * @param n the direction vector
         * @param mode 0 for gluon, +/-1 for quark, +/-2 for anti-quark, + for out, - form in
         * @return Eikonal-Gloun Vertex
         */
        ex eikonalVertex(ex e, ex n, int mode) {
            auto fi1 = e.op(0).op(1);
            auto fi2 = e.op(1).op(1);
            auto fi3 = e.op(2).op(1);
            auto mom1 = e.op(0).op(2);
            if(mode==0) return I * gs * SP(n,LI(fi3)) * (-I*SUNF(CI(fi3),CI(fi1),CI(fi2)));
            else if(mode==1 || mode==-2) return I * gs * SP(n,LI(fi3)) * GMat(GAS(1), DI(fi1),DI(fi2)) * (-SUNT(CI(fi3),TI(fi2),TI(fi1)));
            else if(mode==2 || mode==-1) return I * gs * SP(n,LI(fi3)) * GMat(GAS(1), DI(fi1),DI(fi2)) * SUNT(CI(fi3),TI(fi1),TI(fi2));
            else return 0;
        }
        
        /**
         * @brief Eikonal-Gloun Vertex, right side w.r.t. cutting line
         * @param e expression with head of Vertex
         * @param n the direction vector
         * @param mode 0 for gluon, +/-1 for quark, +/-2 for anti-quark, + for out, - form in
         * @return Eikonal-Gloun Vertex
         */
        ex eikonalVertexR(ex e, ex n, int mode) {
            auto fi1 = e.op(0).op(1);
            auto fi2 = e.op(1).op(1);
            auto fi3 = e.op(2).op(1);
            auto mom1 = e.op(0).op(2);
            if(mode==0) return -I * gs * SP(n,LI(fi3)) * GMat(GAS(1), DI(fi1),DI(fi2)) * (I*SUNF(CI(fi3),CI(fi1),CI(fi2)));
            else if(mode==1 || mode==-2) return -I * gs * SP(n,LI(fi3)) * GMat(GAS(1), DI(fi1),DI(fi2)) * (-SUNT(CI(fi3),TI(fi2),TI(fi1)));
            else if(mode==2 || mode==-1) return -I * gs * SP(n,LI(fi3)) * GMat(GAS(1), DI(fi1),DI(fi2)) * SUNT(CI(fi3),TI(fi1),TI(fi2));
            else return 0;
        }
        
        /**
         * @brief Gluon Fragmentation Function Vertex, nbar-e-g
         * @param e expression with head of Vertex
         * @param n the direction vector
         * @return Eikonal-Gloun Vertex
         */
         ex GluonFFV(ex e, ex n) {
            auto fi1 = e.op(0).op(1);
            auto fi2 = e.op(1).op(1);
            auto fi3 = e.op(2).op(1);
            auto mom2 = e.op(1).op(2);
            auto mom3 = e.op(2).op(2);
            return I * (SP(n,mom2)*SP(LI(fi2),LI(fi3)) + SP(mom3,LI(fi2))*SP(n,LI(fi3))) * SP(CI(fi1),CI(fi3));
         }
              
       /**
        * @brief Quark Fragmentation Function Vertex, qbar-e-nbar/Qbar-e-nbar
        * @param e expression with head of Vertex
        * @param n the direction vector
        * @return Eikonal-Gloun Vertex
        */
        ex QuarkFFV(ex e, ex n) {
            auto fi1 = e.op(0).op(1);
            auto fi3 = e.op(2).op(1);
            return SP(TI(fi1),TI(fi3)) * GMat(GAS(1), DI(fi1),DI(fi3));
        }
    
        // mode = 0 for gluon
        lst FeynRules(const lst & amps, const ex & m0, int mode) {
            lst ret;
            for(auto item : amps) ret.append(FeynRules(item,m0,mode));
            return ret;
        }
    
        // mode = 0 for gluon
        ex FeynRules(const ex & amp, const ex & m0, int mode) {
            if(is_a<lst>(amp)) return FeynRules(ex_to<lst>(amp),m0,mode);
            static Symbol ep("e"), n("n"), nbar("nbar");
            static Symbol Q("Q"), Qbar("Qbar"), q("q"), qbar("qbar"), g("g"), gh("gh"), ghbar("ghbar");
            static Vector vn("n"), q1("q1"), q2("q2"), q3("q3");
            auto fr = MapFunction([mode,m0](const ex &e, MapFunction &self)->ex {
                if(isFunction(e,"OutField") || isFunction(e,"InField")) return 1;
                else if(isFunction(e, "Propagator")) {
                    if(e.op(0).op(0)==q) {
                        return QuarkPropagator(e, 0);
                    } else if(e.op(0).op(0)==Q) {
                        return QuarkPropagator(e, m0);
                    } else if(e.op(0).op(0)==g) {
                        return GluonPropagator(e);
                    } else if(e.op(0).op(0)==gh) {
                        return GluonGhostPropagator(e);
                    } else if(e.op(0).op(0)==n) {
                        auto ret = eikonalPropagator(e, vn, mode);
                        if(!ret.has(q1) && !ret.has(q2) && !ret.has(q3)) ret = ret.subs(iEpsilon==0);
                        return ret;
                    }
                } else if(isFunction(e, "Vertex")) {
                    if(e.nops()==3 && e.op(0).op(0)==nbar && e.op(1).op(0)==ep && e.op(2).op(0)==g) {
                        // nbar-e-g
                        return GluonFFV(e, vn);
                    } else if(e.nops()==3 && e.op(0).op(0)==Qbar && e.op(1).op(0)==ep && e.op(2).op(0)==nbar) {
                        // Qbar-e-nbar
                        return QuarkFFV(e, vn);
                    } else if(e.nops()==3 && e.op(0).op(0)==nbar && e.op(1).op(0)==n && e.op(2).op(0)==g) {
                        // nbar-n-g
                        return eikonalVertex(e, vn, mode);
                    } else if(e.nops()==3 && ((e.op(0).op(0)==qbar && e.op(1).op(0)==q) || (e.op(0).op(0)==Qbar && e.op(1).op(0)==Q)) && (e.op(2).op(0)==g) ) {
                        // qbar-q-g or Qbar-Q-g
                        return QuarkGluonVertex(e);
                    } else if(e.nops()==3 && e.op(0).op(0)==ghbar && e.op(1).op(0)==gh) {
                        // ghbar-gh-g
                        return GhostGluonVertex(e);
                    } else if(e.nops()==3 && e.op(0).op(0)==g && e.op(1).op(0)==g && e.op(2).op(0)==g) {
                        // g^3
                        return Gluon3Vertex(e);
                    } else if(e.nops()==4 && e.op(0).op(0)==g && e.op(1).op(0)==g && e.op(2).op(0)==g && e.op(3).op(0)==g) {
                        // g^4
                        return Gluon4Vertex(e);
                    }
                } else return e.map(self);
                return e;
            });
            return fr(amp);
        }
        
        // mode = 0 for gluon
        ex eSUM(int mode) {
            if(mode==0) return SP(LI(-1), RLI(-1)) * SP(CI(-2),RCI(-2));
            else return SP(TI(-2),RTI(-2)) * GMat(GAS(Vector("n")), DI(-2), RDI(-2));
        }
        
        /**
         * @brief the prefactor for the z-Integration
         * arXiv:1208.5301v2 [hep-ph] 9 Nov 2012
         * PHYSICAL REVIEW D 91, 074013 (2015)
         * @param tls the number of particles in cut
         * @param SF the symmetry factor n, will used as 1/n
         * @param mode 0-Gluon, 1-Quark
         * @return the prefactor for the z-Integration
         */
        ex zIntFactor(int tls, const ex SF, int mode) {
            Symbol kp("kp"), z("z"), M("M");
            ex ret = 0, NCS = 1;
            if(mode==0) NCS = -pow(z,1-2*ep)/(16*(2-2*ep)*kp*Pi);
            else if(mode==1) NCS = pow(z,(1-2*ep))/(8*3*Pi);
            else throw Error("zIntFactor: Not supported mode "+to_string(mode));
            
            if(tls>0) ret = NCS*M/SF*(-((pow(2,(2-tls-(3-2*ep)*tls))*pow(Pi,(1-(3-2*ep)*tls)))/(kp*(-1+z))));
            else ret = NCS*M*4*Pi/kp;
            return Symbol::set_all(ret);
        }
        
        /**
         * @brief dPhi (c1 p.k + c0)^(-n), NOT including the zIntFactor yet
         * @param c1 coefficient of p.k
         * @param c0 constant term
         * @param n the exponent
         * @param k2 the k^2 = k.k
         * @param p2 the p^2 = p.p
         * @return the z integration
         */
        ex zIntegrate(const ex & c1, const ex & c0, const ex & n, const ex & k2, const ex& p2) {
            Symbol z("z"), kp("kp"), pp("pp"); // pp = p.n, 1-z = k.n/kp
            ex ret =  pow(2,n)*pow(Pi,1-ep)*pow(2*c0+c1*kp*(p2-p2*z)*pow(pp,-1)+c1*k2*pp*pow(kp-kp*z,-1),-n)*pow(c1*pow(pp,2)*pow(-2*c0*kp*pp*(-1+z)+c1*(k2*pow(pp,2)+p2*pow(kp,2)*pow(-1+z,2)),-1),-1+ep)*pow(tgamma(n),-1)*tgamma(-1+ep+n);
            return Symbol::set_all(ret);
        }
        
    }
}

