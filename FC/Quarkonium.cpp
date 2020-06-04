/**
 * @file
 * @brief Helpers for Quarkonium
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-20
 */
 
#include "FC.h"
#include "SD.h"

namespace HepLib::FC {
    using namespace Qgraf;

    //-----------------------------------------------------------
    // Quarkonium Class
    //-----------------------------------------------------------
    namespace Quarkonium {
    
        ex Gamma5(const string pre, int start) {
            Index i1(pre+to_string(start+0), Index::Type::VD);
            Index i2(pre+to_string(start+1), Index::Type::VD);
            Index i3(pre+to_string(start+2), Index::Type::VD);
            Index i4(pre+to_string(start+3), Index::Type::VD);
            return LC(i1,i2,i3,i4)*GAS(i1)*GAS(i2)*GAS(i3)*GAS(i4)/(I*factorial(4));
        }
            
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
            if(io==IO::Out) return Matrix(SpinProj(io,s,p,pb,m,e,mu), Qgraf::DI(j), Qgraf::DI(i));
            else return Matrix(SpinProj(io,s,p,pb,m,e,mu), Qgraf::DI(i), Qgraf::DI(j));
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
            return Matrix(SpinProj(io,s,p,pb,m,e,mb,eb,mu), Qgraf::DI(j), Qgraf::DI(i));
        }
        
        /**
         * @brief Color-Singlet Projector 
         * @param i quark qgraf index
         * @param j anti-quark qgraf index
         * @return Color-Singlet Projector
         */
        ex ColorProj(int i, int j) {
            return SP(Qgraf::TI(i),Qgraf::TI(j))/sqrt(NF);
        }
        
        /**
         * @brief Color-Octet Projector 
         * @param i quark qgraf index
         * @param j anti-quark qgraf index
         * @return Color-Octet Projector
         */
        ex ColorProj(int i, int j, Index a) {
            return sqrt(ex(2)) * SUNT(Qgraf::TI(i), Qgraf::TI(j), a);
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
        ex SL1Proj(ex si, ex qi, ex p) {
            return ITD(si,qi,p)/sqrt(D-1);
        }
        
        /**
         * @brief S-L with total spin 1, L=1
         * @param si spin index
         * @param qi q-vector index
         * @param mu spin index
         * @param total momentum
         * @return S-L with total spin 1 
         */
        ex SL1Proj(ex si, ex qi, ex mu, ex p) {
            return -I*LC(si,qi,mu,p)/sqrt(2*SP(p));
        }
        
        /**
         * @brief S-L with total spin 2, L=1
         * @param si spin index
         * @param qi q-vector index
         * @param mu1 spin-2 index
         * @param mu2 spin-2 index
         * @param total momentum
         * @return S-L with total spin 2
         */
        ex SL1Proj(ex si, ex qi, ex mu1, ex mu2, ex p) {
            return (ITD(si,mu1,p)*ITD(qi,mu2,p)+ITD(qi,mu1,p)*ITD(si,mu2,p))/2 - ITD(si,qi,p)*ITD(mu1,mu2,p)/(D-1);
        }
        
        /**
         * @brief S-L with total spin 1, L=2
         * @param si spin index
         * @param qi1 q1-vector index
         * @param qi2 q2-vector index
         * @param mu spin index
         * @param total momentum
         * @return S-L with total spin 1
         */
        ex SL2Proj(ex si, ex qi1, ex qi2, ex mu, ex p) {
            return sqrt((D-1)/(D+1))*((ITD(si,qi1,p)*ITD(qi2,mu,p)+ITD(si,qi2,p)*ITD(qi1,mu,p))/2- ITD(si,mu,p)*ITD(qi1,qi2,p)/(D-1));
        }
        
        /**
         * @brief S-L with total spin 2, L=2
         * @param si spin index
         * @param qi1 q1-vector index
         * @param qi2 q2-vector index
         * @param mu1 spin-2 index
         * @param mu2 spin-2 index
         * @param total momentum
         * @return S-L with total spin L
         */
        ex SL2Proj(ex si, ex qi1, ex qi2, ex mu1, ex mu2, ex p) {
            return -I/(sqrt(2*(D-1)*SP(p))) * (SP(mu1,qi1)*LC(qi2,si,mu2,p)+SP(mu1,qi2)*LC(qi1,si,mu2,p));
        }
        
        /**
         * @brief S-L spin sum for L=0,1,2
         * @param si spin index
         * @param siR spin index @ right side
         * @param qi q-vector index
         * @param qiR q-vector index @ right side
         * @param p total momentum
         * @param L total spin, L=0,1,2
         * @return S-L with total spin L
         */
        ex SLSum(ex si, ex siR, ex qi, ex qiR, ex p, int L) {
            if(L==0) return ITD(si,qi,p)*ITD(siR,qiR,p)/(D-1);
            else if(L==1) return (ITD(si,siR,p)*ITD(qi,qiR,p)-ITD(si,qiR,p)*ITD(qi,siR,p))/2;
            else if(L==2) return (ITD(si,siR,p)*ITD(qi,qiR,p)+ITD(si,qiR,p)*ITD(qi,siR,p))/2-ITD(si,qi,p)*ITD(siR,qiR,p)/(D-1);
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
            expr = mma_collect(expr, q);
            int lproj=0;

            expr = MapFunction([&lproj,prefix,q](const ex &e, MapFunction &self)->ex {
                if(!e.has(q)) return e;
                else if(is_a<Pair>(e)) {
                    if(is_a<Index>(e.op(1))) return e;
                    Index idx(prefix+to_string(lproj++));
                    return SP(e.op(0), idx) * SP(e.op(1), idx);
                } else if(is_a<Eps>(e)) {
                    Index idx(prefix+to_string(lproj++));
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
                    Index idx(prefix+to_string(lproj++));
                    auto g = ex_to<DiracGamma>(e);
                    return DiracGamma(idx, g.rl) * SP(g.pi, idx);
                } else if (e.match(TR(w))) {
                    auto ret = self(e.op(0));
                    ret = mma_collect(ret, q, true);
                    ret = ret.subs(coCF(w)==TR(w));
                    return ret;
                } else if(is_a<add>(e)) {
                    auto lpj = lproj;
                    ex ret = 0;
                    for(auto item : e) {
                        lproj = lpj;
                        ret += self(item);
                    }
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
            auto cv_lst = mma_collect_lst(expr, q);
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
                        is.push_back(ex_to<Index>(item.op(1)));
                    } else if(item.match(pow(w,2)) && is_a<Pair>(item.op(0))) {
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
                        else if(isn==2) qproj = SP(q)*(SP(is[0],is[1])-SP(p,is[0])*SP(p,is[1])/SP(p))/(D-1);
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
                            qproj = -SP(q)*(ITD(m1,m2,p)*SP(qi,m3)+ITD(m1,m3,p)*SP(qi,m2)+ITD(m2,m3,p)*SP(qi,m1))/(D+1);
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
        
        /**
         * @brief n-massless body phase space
         * https://arxiv.org/abs/hep-ph/0311276v1
         * @param n the number of massless particles
         * @return integrated phase space
         */
        ex nPS(int n, ex q2) {
            return pow(2,5-4*n-2*ep+2*n*ep) * pow(Pi,3-2*n-ep+n*ep) *
                tgamma(1-ep)/(tgamma((n-1)*(1-ep))*tgamma(n*(1-ep))) * pow(q2,n-2+ep-n*ep);
        }
        
        namespace {
            inline ex V(ex dim) { return 2*pow(Pi,dim/2) / tgamma(dim/2); }
        }
        
        /**
         * @brief 2-, 3-, 4- massless Phase space
         * https://arxiv.org/abs/hep-ph/0311276v1
         * @param moms momentum in phase space
         * @param amp input amplitudes
         * @param si start index, -1 for automatically detect
         * @param q2 refers to q^2, total invarant mass 
         * @return amp * PS2/3/4, for PS4, a list will be returned
         */
        ex DoPS(lst moms, ex amp, int si_in, ex q2) {
            int si = si_in;
            if(si<0) {
                auto xs = SD::get_x_from(amp);
                if(xs.size()<1) si=0;
                else si=ex_to<numeric>(xs[xs.size()-1].op(0)).to_int()+1;
            }
            auto d = 4-2*ep;
            int nc = moms.nops();
            if(nc==2) {
                return pow(2*Pi,2-d) * pow(q2,(d-4)/2) * V(d-1)/pow(2,d-1) * amp;
            } else if(nc==3) {
                auto p1 = moms.op(0);
                auto p2 = moms.op(1);
                auto p3 = moms.op(2);
                ex s12 = x(si+0)*q2;
                ex s13 = x(si+1)*q2;
                ex s23 = x(si+2)*q2;
                exmap sp2x;
                sp2x[p1*p2]=s12/2;
                sp2x[p1*p3]=s13/2;
                sp2x[p2*p3]=s23/2;
                sp2x[p1*p1]=0;
                sp2x[p2*p2]=0;
                sp2x[p3*p3]=0;

                auto ret = pow(2*Pi,3-2*d) * pow(2,-1-d) * pow(q2,(2-d)/2) * V(d-1)*V(d-2) * amp.subs(sp2x,subs_options::algebraic);

                // add (s12 s13 s23)^((d-4)/2) to each F
                ex ss = s12*s13*s23;
                ret = MapFunction([si,ss,d](const ex & e, MapFunction &self)->ex{
                    if(!XIntegral::has(e)) return e;
                    else if(is_a<XIntegral>(e)) {
                        XIntegral xint = ex_to<XIntegral>(e);
                        let_op_append(xint.Functions,ss);
                        let_op_append(xint.Exponents,(d-4)/2);
                        let_op_append(xint.Deltas,lst{x(si+0),x(si+1),x(si+2)});
                        return xint;
                    } else return e.map(self);
                })(ret);
                
                return ret;
            } else if(nc==4) {
                auto p1 = moms.op(0);
                auto p2 = moms.op(1);
                auto p3 = moms.op(2);
                auto p4 = moms.op(3);
                
                ex s12 = x(si+0)*q2;
                ex s13 = x(si+1)*q2;
                ex s23 = x(si+2)*q2;
                ex s14 = x(si+3)*q2;
                ex s24 = x(si+4)*q2;
                ex s34 = x(si+5)*q2;
                exmap sp2x;
                sp2x[p1*p1]=0;
                sp2x[p2*p2]=0;
                sp2x[p3*p3]=0;
                sp2x[p4*p4]=0;
                sp2x[p1*p2]=s12/2;
                sp2x[p1*p3]=s13/2;
                sp2x[p2*p3]=s23/2;
                sp2x[p1*p4]=s14/2;
                sp2x[p2*p4]=s24/2;
                sp2x[p3*p4]=s34/2;
                
                auto ret = pow(2*Pi,4-3*d) * pow(2,1-2*d) * V(d-1)*V(d-2)*V(d-3) *
                    pow(q2,3*d/2-4) * amp.subs(sp2x,subs_options::algebraic);
                
                // add λ(x16,x25,x34)^((d-5)/2) to each F
                ret = MapFunction([si,d](const ex & e, MapFunction &self)->ex{
                    if(!XIntegral::has(e)) return e;
                    else if(is_a<XIntegral>(e)) {
                        XIntegral xint = ex_to<XIntegral>(e);
                        ex lambda = pow(x(si+2),2)*pow(x(si+3),2) + pow(x(si+1)*x(si+4)-x(si+0)*x(si+5),2) - 2*x(si+2)*x(si+3)*(x(si+1)*x(si+4) + x(si+0)*x(si+5));
                        let_op_append(xint.Functions,-lambda);
                        let_op_append(xint.Exponents,(d-5)/2); 
                                                
                        lst xs;
                        for(int i=0; i<6; i++) xs.append(x(si+i));
                        let_op_append(xint.Deltas,xs);
                        
                        ex fe = lst{xint.Functions, xint.Exponents, xint.Deltas};
                        SD::SecDec::Projectivize(fe,xs);
                        
                        // cheng-wu: change x(0,1,2) from 0 to ∞, and rescale
                        SD::SecDec::Scalelize(fe,x(si+0),1/x(si+5));
                        SD::SecDec::Scalelize(fe,x(si+1),1/x(si+4));
                        SD::SecDec::Scalelize(fe,x(si+2),1/x(si+3));
                        auto fe_vec = SD::SecDec::Binarize(fe, x(si+0)-x(si+2));
                        auto fe0 = fe_vec[0];
                        auto fe2 = fe_vec[1];
                        
                        auto fe0_vec = SD::SecDec::Binarize(fe0, x(si+0)-x(si+1));
                        auto fe2_vec = SD::SecDec::Binarize(fe2, x(si+2)-x(si+1));
                        
                        SD::SecDec::Scalelize(fe0_vec[0],x(si+2),pow(x(si+0)/x(si+1),2));
                        // now Θ refers to x(2)>x(1)
                        fe0_vec[0] = SD::SecDec::Binarize(fe0_vec[0], x(si+2)-x(si+1))[0];
                        
                        SD::SecDec::Scalelize(fe0_vec[1],x(si+2),pow(x(si+1),2)/(x(si+0)+2*x(si+1))/x(si+0));
                        // now Θ refers to x(2)>x(0)
                        fe0_vec[1] = SD::SecDec::Binarize(fe0_vec[1], x(si+2)-x(si+0))[0];

                        SD::SecDec::Scalelize(fe2_vec[0],x(si+0),pow(x(si+2)/x(si+1),2));
                        // now Θ refers to x(0)>x(1)
                        fe2_vec[0] = SD::SecDec::Binarize(fe2_vec[0], x(si+0)-x(si+1))[0];
                        
                        SD::SecDec::Scalelize(fe2_vec[1],x(si+0),pow(x(si+1),2)/(x(si+2)+2*x(si+1))/x(si+2));
                        // now Θ refers to x(0)>x(2)
                        fe2_vec[1] = SD::SecDec::Binarize(fe2_vec[1], x(si+0)-x(si+2))[0];
                        
                        ex xsum = x(si+3)+x(si+4)+x(si+5);
                        xsum = x(si+5);
                        SD::SecDec::Projectivize(fe0_vec[0],xs,xsum);
                        SD::SecDec::Projectivize(fe0_vec[1],xs,xsum);
                        SD::SecDec::Projectivize(fe2_vec[0],xs,xsum);
                        SD::SecDec::Projectivize(fe2_vec[1],xs,xsum);
                        
                        return XIntegral(fe0_vec[0]) + XIntegral(fe0_vec[1]) + XIntegral(fe2_vec[0]) + XIntegral(fe2_vec[1]);
                    } else return e.map(self);
                })(ret);
                
                return ret;
            } else throw Error("DoPS: only 2-, 3-, 4- massless Phase Space supported.");
        }
        
    }
}

