/**
 * @file
 * @brief Helpers for Quarkonium
 */

#include "QCD.h"
#include "SD.h"

namespace HepLib::QCD {

    //-----------------------------------------------------------
    // Quarkonium Class
    //-----------------------------------------------------------
    namespace Quarkonium {
    
        using namespace SD;
            
        /**
         * @brief Spin Projector with the same mass, non-Relativistic Normalization, 
         * To Use Relativistic Normalization, Just Times 2e
         * SpinProj[In,p,pb,m,e,mu], where p=P/2+q, pb=P/2-q, e^2=m^2-q^2,
         * In CM, P=(M=2e,0) & q=(0,q), P.q=0, P.q=0
         * @param io In or Out
         * @param s s=0 or s=1
         * @param p quark momentum
         * @param pb anti-quark momentum
         * @param m mass
         * @param e e^2=m^q-q^2
         * @param mu only for s=1
         * @return the corresponding spin projector
         */
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mu) {
            ex ret;
            if(io==In && s==1) ret = (GAS(p)+m*GAS(1))*(GAS(p+pb)+2*e*GAS(1))*GAS(mu)*(GAS(pb)-m*GAS(1));
            if(io==Out && s==1) ret = (GAS(pb)-m*GAS(1))*GAS(mu)*(GAS(p+pb)+2*e*GAS(1))*(GAS(p)+m*GAS(1));
            if(io==In && s==0) ret = -(GAS(p)+m*GAS(1))*(GAS(p+pb)+2*e*GAS(1))*GAS(5)*(GAS(pb)-m*GAS(1));
            if(io==Out && s==0) ret = (GAS(pb)-m*GAS(1))*GAS(5)*(GAS(p+pb)+2*e*GAS(1))*(GAS(p)+m*GAS(1));
            ret = ret/(4*sqrt(ex(2))*e*(e+m));
            ret = ret/(2*e);
            return ret;
        }
        
        /**
         * @brief Spin Projector with the same mass, non-Relativistic Normalization, 
         * - To Use Relativistic Normalization, Just Times 2e
         * - SpinProj[In,p,pb,m,e,mu], where p=P/2+q, pb=P/2-q, e^2=m^2-q^2,
         * - In CM, P=(M=2e,0) & q=(0,q), P.q=0, P.q=0
         * @param io In or Out
         * @param s s=0 or s=1
         * @param p quark momentum
         * @param pb anti-quark momentum
         * @param m mass
         * @param e e^2=m^q-q^2
         * @param mu only for s=1
         * @param i quark qgraf index
         * @param j anti-quark qgraf index
         * @return the corresponding spin projector
         */
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mu, int i, int j) {
            if(io==IO::Out) return Matrix(SpinProj(io,s,p,pb,m,e,mu), DI(j), DI(i));
            else return Matrix(SpinProj(io,s,p,pb,m,e,mu), DI(i), DI(j));
        }
        
        /**
         * @brief Spin Projector with different mass, non-Relativistic Normalization, 
         * - To Use Relativistic Normalization, Just Times Sqrt[2e 2eb]
         * - SpinProj[In,p,pb,m,mb,e,eb,mu], where p=e/(e+eb)P+q,pb=eb/(e+eb)P-q,e^2=m^2-q^2,eb^2=mb^2-qb^2,
         * - In CM, P=(M=e+eb,0) & q=(0,q), P.q=0
         * - arXiv:1003.0061
         * @param io In or Out
         * @param s s=0 or s=1
         * @param p quark momentum
         * @param pb anti-quark momentum
         * @param m quark mass
         * @param mb anti-quark mass
         * @param e for quark: e^2=m^q-q^2
         * @param eb for anti-quark eb^2=mb^2-qb^2 
         * @param mu only for s=1
         * @return the corresponding spin projector
         */
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mb, ex eb, ex mu) {
            ex ret;
            if(io==In && s==1) ret = (GAS(p)+m*GAS(1))*(GAS(p+pb)+(e+eb)*GAS(1))*GAS(mu)*(GAS(pb)-mb*GAS(1));
            if(io==Out && s==1) ret = (GAS(pb)-mb*GAS(1))*GAS(mu)*(GAS(p+pb)+(e+eb)*GAS(1))*(GAS(p)+m*GAS(1));
            if(io==In && s==0) ret = -(GAS(p)+m*GAS(1))*(GAS(p+pb)+(e+eb)*GAS(1))*GAS(5)*(GAS(pb)-mb*GAS(1));
            if(io==Out && s==0) ret = (GAS(pb)-mb*GAS(1))*GAS(5)*(GAS(p+pb)+(e+eb)*GAS(1))*(GAS(p)+m*GAS(1));
            ret = ret/(2*sqrt(ex(2))*(e+eb)*sqrt((e+m)*(eb+mb)));
            ret = ret/(sqrt(4*e*eb));
            return ret;
        }
        
        /**
         * @brief Spin Projector with different mass, non-Relativistic Normalization, 
         * - To Use Relativistic Normalization, Just Times Sqrt[2e 2eb]
         * - SpinProj[In,p,pb,m,mb,e,eb,mu], where p=e/(e+eb)P+q,pb=eb/(e+eb)P-q,e^2=m^2-q^2,eb^2=mb^2-qb^2,
         * - In CM, P=(M=e+eb,0) & q=(0,q), P.q=0
         * - arXiv:1003.0061
         * @param io In or Out
         * @param s s=0 or s=1
         * @param p quark momentum
         * @param pb anti-quark momentum
         * @param m quark mass
         * @param mb anti-quark mass
         * @param e for quark: e^2=m^q-q^2
         * @param eb for anti-quark eb^2=mb^2-qb^2 
         * @param mu only for s=1
         * @param i quark qgraf index
         * @param j anti-quark qgraf index
         * @return the corresponding spin projector
         */
        ex SpinProj(IO io, int s, ex p, ex pb, ex m, ex e, ex mb, ex eb, ex mu, int i, int j) {
            return Matrix(SpinProj(io,s,p,pb,m,e,mb,eb,mu), DI(j), DI(i));
        }
        
        /**
         * @brief Color-Singlet Projector 
         * @param i quark qgraf index
         * @param j anti-quark qgraf index
         * @return Color-Singlet Projector
         */
        ex ColorProj(int i, int j) {
            return SP(TI(i),TI(j))/sqrt(NF);
        }
        
        /**
         * @brief Color-Octet Projector 
         * @param i quark qgraf index
         * @param j anti-quark qgraf index
         * @param a the color index
         * @return Color-Octet Projector
         */
        ex ColorProj(int i, int j, Index a) {
            return sqrt(ex(2)) * SUNT(a, TI(i), TI(j));
        }
        
        namespace {
            ex ITD(ex si, ex qi , ex p) {
                return -SP(si,qi)+SP(p,si)*SP(p,qi)/SP(p);
            }
        }
        
        /**
         * @brief S-L with total spin 0, L=1
         * @param si spin index
         * @param qi q-vector index
         * @param p total momentum
         * @return S-L with total spin 0
         */
        ex S1L1Proj(ex si, ex qi, ex p) {
            return ITD(si,qi,p)/sqrt(d-1);
        }
        
        /**
         * @brief S-L with total spin 1, L=1
         * @param si spin index
         * @param qi q-vector index
         * @param mu spin index
         * @param p total momentum
         * @return S-L with total spin 1 
         */
        ex S1L1Proj(ex si, ex qi, ex mu, ex p) {
            return -I*LC(si,qi,mu,p)/sqrt(2*SP(p));
        }
        
        /**
         * @brief S-L with total spin 2, L=1
         * @param si spin index
         * @param qi q-vector index
         * @param mu1 spin-2 index
         * @param mu2 spin-2 index
         * @param p total momentum
         * @return S-L with total spin 2
         */
        ex S1L1Proj(ex si, ex qi, ex mu1, ex mu2, ex p) {
            return (ITD(si,mu1,p)*ITD(qi,mu2,p)+ITD(qi,mu1,p)*ITD(si,mu2,p))/2 - ITD(si,qi,p)*ITD(mu1,mu2,p)/(d-1);
        }
        
        /**
         * @brief S-L with total spin 1, L=2
         * @param si spin index
         * @param qi1 q1-vector index
         * @param qi2 q2-vector index
         * @param mu spin index
         * @param p total momentum
         * @return S-L with total spin 1
         */
        ex S1L2Proj(ex si, ex qi1, ex qi2, ex mu, ex p) {
            return sqrt((d-1)/(d+1))*((ITD(si,qi1,p)*ITD(qi2,mu,p)+ITD(si,qi2,p)*ITD(qi1,mu,p))/2- ITD(si,mu,p)*ITD(qi1,qi2,p)/(d-1));
        }
        
        /**
         * @brief S-L with total spin 2, L=2
         * @param si spin index
         * @param qi1 q1-vector index
         * @param qi2 q2-vector index
         * @param mu1 spin-2 index
         * @param mu2 spin-2 index
         * @param p total momentum
         * @return S-L with total spin L
         */
        ex S1L2Proj(ex si, ex qi1, ex qi2, ex mu1, ex mu2, ex p) {
            return -I/(sqrt(2*(d-1)*SP(p))) * (SP(mu1,qi1)*LC(qi2,si,mu2,p)+SP(mu1,qi2)*LC(qi1,si,mu2,p));
        }
        
        /**
         * @brief S-L spin sum for J=0,1,2
         * @param si spin index
         * @param siR spin index @ right side
         * @param qi q-vector index
         * @param qiR q-vector index @ right side
         * @param p total momentum
         * @param J total spin, J=0,1,2
         * @return S-L with total spin J
         */
        ex S1L1Sum(ex si, ex siR, ex qi, ex qiR, ex p, int J) {
            if(J==0) return ITD(si,qi,p)*ITD(siR,qiR,p)/(d-1);
            else if(J==1) return (ITD(si,siR,p)*ITD(qi,qiR,p)-ITD(si,qiR,p)*ITD(qi,siR,p))/2;
            else if(J==2) return (ITD(si,siR,p)*ITD(qi,qiR,p)+ITD(si,qiR,p)*ITD(qi,siR,p))/2-ITD(si,qi,p)*ITD(siR,qiR,p)/(d-1);
            return 0;
        }
        
        /**
         * @brief S/P/D Wave Project
         * @param expr_in the input expression
         * @param pqi the lst { p, q, i1,i2,i3,... }
         * @param prefix the prefix used in the uncontract dumy index
         * @return the LProj result
         */
        ex LProj(const ex &expr_in, const lst &pqi, string prefix) {
            ex p = pqi.op(0);
            ex q = pqi.op(1);
            
            if(!is_a<Vector>(q)) throw Error("LProj invalid 3rd argument");
            if(expr_in.has(coVF(w))) throw Error("LProj error: expr_in has coVF already.");
            for(int i=2; i<pqi.nops(); i++) {
                if(!is_a<Index>(pqi.op(i))) throw Error("LProj invalid 3rd argument");
            }
            
            // Un-Contract
            auto expr = expr_in.subs(SP_map);
            expr = collect_ex(expr, q);
            int lproj=-1;

            expr = MapFunction([&lproj,prefix,q](const ex &e, MapFunction &self)->ex {
                if(!e.has(q)) return e;
                else if(is_a<Pair>(e)) {
                    if(!e.op(0).is_equal(q) && !e.op(1).is_equal(q)) throw Error("LProj: e.op(0) is NOT q with e = " + ex2str(e));
                    if(is_a<Index>(e.op(1))) return e;
                    Index idx(prefix+to_string(++lproj));
                    return SP(e.op(0), idx) * SP(e.op(1), idx);
                } else if(is_a<Eps>(e)) {
                    auto pis0 = ex_to<Eps>(e).pis;
                    ex pis[4];
                    ex cc = 1;
                    for(int i=0; i<4; i++) {
                        pis[i] = pis0[i];
                        if(is_zero(pis[i]-q)) {
                            Index idx(prefix+to_string(++lproj));
                            pis[i] = idx;
                            cc *= SP(q, idx);
                        } else if(pis[i].has(q)) throw Error("LProj: Eps still has q.");
                    }
                    return LC(pis[0], pis[1], pis[2], pis[3]) * cc;
                } else if(is_a<DGamma>(e)) {
                    Index idx(prefix+to_string(++lproj));
                    auto g = ex_to<DGamma>(e);
                    if(!g.pi.is_equal(q)) throw Error("LProj: g.pi is NOT q.");
                    return DGamma(idx, g.rl) * SP(g.pi, idx);
                } else if (e.match(TR(w))) {
                    auto ret = self(e.op(0));
                    ret = collect_ex(ret, q, true);
                    ret = ret.subs(coCF(w)==TR(w));
                    return ret;
                } else if(is_a<add>(e)) {
                    int lpj = lproj;
                    int lpj_max = -100;
                    ex ret = 0;
                    for(auto item : e) {
                        lproj = lpj;
                        ret += self(item);
                        if(lpj_max<lproj) lpj_max = lproj;
                    }
                    lproj = lpj_max;
                    return ret;
                } else if(is_a<power>(e)) {
                    if(!e.op(1).info(info_flags::posint)) {
                        cout << e << endl;
                        throw Error("LProj: power is not info_flags::posint.");
                    }
                    ex ret = 1;
                    int pn = ex_to<numeric>(e.op(1)).to_int();
                    for(int i=0; i<pn; i++) {
                        ret *= self(e.op(0));
                    }
                    return ret;
                } else return e.map(self);
                throw Error("something should be wrong here");
                return 0;
            })(expr);
            
            expr = expr.subs(SP_map);
            auto cv_lst = collect_lst(expr, q);

            expr = 0;
            for(auto cv : cv_lst) {
                auto e = cv.op(1);
                if(is_zero(e-1)) {
                    if(pqi.nops()==2) expr += cv.op(0);
                    continue;
                }
                
                if(!is_a<mul>(e)) e = lst{e};
                vector<Index> is;
                
                for(auto item : e) {
                    if(is_a<Pair>(item)) {
                        if(!item.op(0).is_equal(q)) throw Error("LProj: op(0) is NOT q.");
                        is.push_back(ex_to<Index>(item.op(1)));
                    } else if(item.match(pow(w,2)) && is_a<Pair>(item.op(0))) {
                        if(!item.op(0).op(0).is_equal(q)) throw Error("LProj: op(0) is NOT q.");
                        is.push_back(ex_to<Index>(item.op(0).op(1)));
                        is.push_back(ex_to<Index>(item.op(0).op(1)));
                    } else {
                        cout << item << endl;
                        throw Error("LProj: something is wrong, unhandled terms.");
                    }
                }

                ex qproj;
                int isn = is.size();
                switch (pqi.nops()) {
                    case 2: {
                        if(isn%2!=0) qproj = 0;
                        else if(isn==2) qproj = -SP(q)/(d-1)*ITD(is[0],is[1],p);
                        else if(isn==4) qproj = pow(SP(q),2)/((d-1)*(d+1)) * (ITD(is[0],is[1],p)*ITD(is[2],is[3],p)+ITD(is[0],is[2],p)*ITD(is[1],is[3],p)+ITD(is[0],is[3],p)*ITD(is[1],is[2],p));
                        else throw Error("LProj not supported yet in S-wave.");
                        break;
                    } case 3: {
                        auto qi = ex_to<Index>(pqi.op(2));
                        if(isn==1) qproj = SP(qi, is[0]);
                        else if(isn==2) qproj = 0;
                        else if(isn==3) {
                            auto m1 = is[0];
                            auto m2 = is[1];
                            auto m3 = is[2];
                            qproj = -SP(q)*(ITD(m1,m2,p)*SP(qi,m3)+ITD(m1,m3,p)*SP(qi,m2)+ITD(m2,m3,p)*SP(qi,m1))/(d+1);
                        } else throw Error("LProj not supported yet in P-wave.");
                        break;
                    } case 4: {
                        auto e1 = ex_to<Index>(pqi.op(2));
                        auto e2 = ex_to<Index>(pqi.op(3));
                        if(isn==1) qproj = 0;
                        if(isn==2) qproj = SP(is[0],e1)*SP(is[1],e2);
                        else throw Error("LProj not supported yet in D-wave.");
                        break;
                    } default: {
                        throw Error("LProj not supported yet with L>2");
                    }
                }
                
                expr += cv.op(0) * qproj;
            }
            
            expr = expr.subs(SP_map);
            return expr;
        }
        
    }
}

