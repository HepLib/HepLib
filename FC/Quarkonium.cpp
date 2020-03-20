#include "FC.h"

namespace HepLib::FC {

    //-----------------------------------------------------------
    // Quarkonium Class
    //-----------------------------------------------------------
    namespace Quarkonium {
    
        // Using non-Relativistic Normalization
        // To Use Relativistic Normalization, Just Times Sqrt[2e 2eb]
        // SProj0In[pc,pcb,m,e,mu], where pc=P/2+q,pcb=P/2-q,e^2=m^2-q^2,
        // In CM, P=(M=2e,0) & q=(0,q), P.q=0, P.q=0
        // SProj0In[pc,pcb,m,mb,e,eb,mu], where p=e/(e+eb)P+q,pb=eb/(e+eb)P-q,e^2=m^2-q^2,eb^2=mb^2-qb^2,
        // In CM, P=(M=e+eb,0) & q=(0,q), P.q=0
        // arXiv:1003.0061
            
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
        
        ex ColorProj() {
            return 1/sqrt(NF);
        }
        ex ColorProj(Index i, Index j, Index a) {
            return sqrt(ex(2)) * SUNT(i,j,a);
        }
        
        namespace {
            ex ITD(ex si, ex qi , ex p) {
                return -SP(si,qi)+SP(p,si)*SP(p,qi)/SP(p);
            }
        }
        
        ex SL1Proj(ex si, ex qi, ex p) {
            return ITD(si,qi,p)/sqrt(D-1);
        }
        
        ex SL1Proj(ex si, ex qi, ex mu, ex p) {
            return -I*LC(si,qi,mu,p)/sqrt(2*SP(p));
        }
        
        ex SL1Proj(ex si, ex qi, ex mu1, ex mu2, ex p) {
            return (ITD(si,mu1,p)*ITD(qi,mu2,p)+ITD(qi,mu1,p)*ITD(si,mu2,p))/2 - ITD(si,qi,p)*ITD(mu1,mu2,p)/(D-1);
        }
        
        ex SL2Proj(ex si, ex qi1, ex qi2, ex mu, ex p) {
            return sqrt((D-1)/(D+1))*((ITD(si,qi1,p)*ITD(qi2,mu,p)+ITD(si,qi2,p)*ITD(qi1,mu,p))/2- ITD(si,mu,p)*ITD(qi1,qi2,p)/(D-1));
        }
        
        ex SL2Proj(ex si, ex qi1, ex qi2, ex mu1, ex mu2, ex p) {
            return -I/(sqrt(2*(D-1)*SP(p))) * (SP(mu1,qi1)*LC(qi2,si,mu2,p)+SP(mu1,qi2)*LC(qi1,si,mu2,p));
        }
        
        ex SLSum(ex si, ex siR, ex qi, ex qiR, ex p, int L) {
            if(L==0) return ITD(si,qi,p)*ITD(siR,qiR,p)/(D-1);
            else if(L==1) return (ITD(si,siR,p)*ITD(qi,qiR,p)-ITD(si,qiR,p)*ITD(qi,siR,p))/2;
            else if(L==2) return (ITD(si,siR,p)*ITD(qi,qiR,p)+ITD(si,qiR,p)*ITD(qi,siR,p))/2-ITD(si,qi,p)*ITD(siR,qiR,p)/(D-1);
            return 0;
        }
        
        // S/P/D Wave Project
        ex LProj(const ex &expr_in, const lst &pqi) {
            static int sidx=0;
            static string prefix = "lproj";
            
            ex p = pqi.op(0);
            ex q = pqi.op(1);
            
            if(!is_a<Vector>(q)) throw Error("LProj invalid 3rd argument");
            if(expr_in.has(coVF(w))) throw Error("LProj error: expr_in has coVF already.");
            for(int i=2; i<pqi.nops(); i++) {
                if(!is_a<Index>(pqi.op(i))) throw Error("LProj invalid 3rd argument");
            }
            
            // Un-Contract
            auto expr = MapFunction([q](const ex &e, MapFunction &self)->ex{
                if(!e.has(q)) return e;
                else if(is_a<Pair>(e)) {
                    if(is_a<Index>(e.op(1))) return e;
                    Index idx(prefix+to_string(sidx++));
                    return SP(e.op(0), idx) * SP(e.op(1), idx);
                } else if(is_a<Eps>(e)) {
                    Index idx(prefix+to_string(sidx++));
                    auto pis0 = ex_to<Eps>(e).pis;
                    ex pis[4];
                    for(int i=0; i<3; i++) {
                        pis[i] = pis0[i];
                        if(is_zero(pis[i]-q)) {
                            pis[i] = idx;
                        }
                    }
                    return LC(pis[0], pis[1], pis[2], pis[3]) * SP(q, idx);
                } else if(is_a<DiracGamma>(e)) {
                    Index idx(prefix+to_string(sidx++));
                    auto g = ex_to<DiracGamma>(e);
                    return DiracGamma(idx, g.rl) * SP(q, idx);
                } else return e.map(self);
                return 0;
            })(expr_in);

            expr = mma_collect(expr, q, false, true);

            expr = MapFunction([q,p,pqi](const ex &e, MapFunction &self)->ex{
                if(e.match(coVF(w))) {
                    vector<Index> is;
                    if(is_a<mul>(e.op(0))) {
                        for(auto item : e.op(0)) is.push_back(ex_to<Index>(item.op(1)));
                    } else if(is_a<Pair>(e.op(0))) {
                        is.push_back(ex_to<Index>(e.op(0).op(1)));
                    } else if(is_zero(e.op(0)-1) && pqi.nops()==2) return 1;
                    else throw Error("something is wrong here.");

                    int isn = is.size();
                    switch (pqi.nops()) {
                        case 2: {
                            if(isn%2!=0) return 0;
                            else if(isn==2) return SP(q)*(SP(is[0],is[1])-SP(p,is[0])*SP(p,is[1])/SP(p))/(D-1);
                            else throw Error("LProj not supported yet in S-wave.");
                            break;
                        } case 3: {
                            auto qi = ex_to<Index>(pqi.op(2));
                            if(isn==1) return SP(qi, is[0]);
                            else if(isn==2) return 0;
                            else if(isn==3) {
                                auto m1 = is[0];
                                auto m2 = is[1];
                                auto m3 = is[2];
                                return -SP(q)*(ITD(m1,m2,p)*SP(qi,m3)+ITD(m1,m3,p)*SP(qi,m2)+ITD(m2,m3,p)*SP(qi,m1))/(D+1);
                            } else throw Error("LProj not supported yet in P-wave.");
                            break;
                        } case 4: {
                            auto e1 = ex_to<Index>(pqi.op(2));
                            auto e2 = ex_to<Index>(pqi.op(3));
                            if(isn==1) return 0;
                            if(isn==2) return SP(is[0],e1)*SP(is[1],e2);
                            else throw Error("LProj not supported yet in D-wave.");
                            break;
                        } default:
                            throw Error("LProj not supported yet with L>2");
                    }
                    return 0; // un-reachable
                } else if (!e.has(coVF(w))) return e;
                else return e.map(self);
            })(expr);
            
            return expr;
        }
        
    }
}

