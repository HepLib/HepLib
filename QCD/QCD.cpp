/**
 * @file
 * @brief Helpers for Quarkonium
 */

#include "QCD.h"
#include "SD.h"

namespace HepLib::QCD {
    
    ex Gamma5(const string pre, int start) {
        Index i1(pre+to_string(start+0), Index::Type::VD);
        Index i2(pre+to_string(start+1), Index::Type::VD);
        Index i3(pre+to_string(start+2), Index::Type::VD);
        Index i4(pre+to_string(start+3), Index::Type::VD);
        return LC(i1,i2,i3,i4)*GAS(i1)*GAS(i2)*GAS(i3)*GAS(i4)/(I*factorial(4));
    }
    
    /**
     * @brief n-massless body phase space
     * https://arxiv.org/abs/hep-ph/0311276v1
     * @param n the number of massless particles
     * @param q2 the squared total momentum
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
     * @param si_in start index, -1 for automatically detect
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
            return pow(2*Pi,2-d) * pow(q2,-ep) * V(d-1)/pow(2,d-1) * amp;
        } else if(nc==3) {
            auto p1 = moms.op(0);
            auto p2 = moms.op(1);
            auto p3 = moms.op(2);
            ex s12 = x(si+0)*q2;
            ex s13 = x(si+1)*q2;
            ex s23 = x(si+2)*q2;
            exmap sp2x;
            sp2x[w*p1*p2]=w*s12/2;
            sp2x[w*p1*p3]=w*s13/2;
            sp2x[w*p2*p3]=w*s23/2;
            sp2x[w*p1*p1]=0;
            sp2x[w*p2*p2]=0;
            sp2x[w*p3*p3]=0;

            auto ret = pow(2*Pi,3-2*d) * pow(2,-1-d) * pow(q2,1-2*ep) * V(d-1)*V(d-2) * amp.subs(sp2x);

            ex ss = x(si+0)*x(si+1)*x(si+2);
            ret = MapFunction([si,ss,d](const ex & e, MapFunction &self)->ex{
                if(!XIntegral::has(e)) return e;
                else if(is_a<XIntegral>(e)) {
                    XIntegral xint = ex_to<XIntegral>(e);
                    let_op_append(xint.Function,ss);
                    let_op_append(xint.Exponent,(d-4)/2);
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
            sp2x[w*p1*p1]=0;
            sp2x[w*p2*p2]=0;
            sp2x[w*p3*p3]=0;
            sp2x[w*p4*p4]=0;
            sp2x[w*p1*p2]=w*s12/2;
            sp2x[w*p1*p3]=w*s13/2;
            sp2x[w*p2*p3]=w*s23/2;
            sp2x[w*p1*p4]=w*s14/2;
            sp2x[w*p2*p4]=w*s24/2;
            sp2x[w*p3*p4]=w*s34/2;
            
            auto ret = pow(2*Pi,4-3*d) * pow(2,1-2*d) * V(d-1)*V(d-2)*V(d-3) *
                pow(q2,2-3*ep) * amp.subs(sp2x);
            
            // add λ(x16,x25,x34)^((d-5)/2) to each F
            ret = MapFunction([si,d](const ex & e, MapFunction &self)->ex{
                if(!XIntegral::has(e)) return e;
                else if(is_a<XIntegral>(e)) {
                    XIntegral xint = ex_to<XIntegral>(e);
                    ex lambda = pow(x(si+2),2)*pow(x(si+3),2) + pow(x(si+1)*x(si+4)-x(si+0)*x(si+5),2) - 2*x(si+2)*x(si+3)*(x(si+1)*x(si+4) + x(si+0)*x(si+5));
                    let_op_append(xint.Function,-lambda);
                    let_op_append(xint.Exponent,(d-5)/2);
                                            
                    lst xs;
                    for(int i=0; i<6; i++) xs.append(x(si+i));
                    let_op_append(xint.Deltas,xs);
                    
                    ex fe = lst{xint.Function, xint.Exponent, xint.Deltas};
                    SD::ChengWu::Projectivize(fe,xs);
                    
                    // cheng-wu: change x(0,1,2) from 0 to ∞, and rescale
                    SD::ChengWu::Scalelize(fe,x(si+0),1/x(si+5));
                    SD::ChengWu::Scalelize(fe,x(si+1),1/x(si+4));
                    SD::ChengWu::Scalelize(fe,x(si+2),1/x(si+3));
                    auto fe_vec = SD::ChengWu::Binarize(fe, x(si+0)-x(si+2));
                    auto fe0 = fe_vec[0];
                    auto fe2 = fe_vec[1];
                    
                    auto fe0_vec = SD::ChengWu::Binarize(fe0, x(si+0)-x(si+1));
                    auto fe2_vec = SD::ChengWu::Binarize(fe2, x(si+2)-x(si+1));
                    
                    SD::ChengWu::Scalelize(fe0_vec[0],x(si+2),pow(x(si+0)/x(si+1),2));
                    // now Θ refers to x(2)>x(1)
                    fe0_vec[0] = SD::ChengWu::Binarize(fe0_vec[0], x(si+2)-x(si+1))[0];
                    
                    SD::ChengWu::Scalelize(fe0_vec[1],x(si+2),pow(x(si+1),2)/(x(si+0)+2*x(si+1))/x(si+0));
                    // now Θ refers to x(2)>x(0)
                    fe0_vec[1] = SD::ChengWu::Binarize(fe0_vec[1], x(si+2)-x(si+0))[0];

                    SD::ChengWu::Scalelize(fe2_vec[0],x(si+0),pow(x(si+2)/x(si+1),2));
                    // now Θ refers to x(0)>x(1)
                    fe2_vec[0] = SD::ChengWu::Binarize(fe2_vec[0], x(si+0)-x(si+1))[0];
                    
                    SD::ChengWu::Scalelize(fe2_vec[1],x(si+0),pow(x(si+1),2)/(x(si+2)+2*x(si+1))/x(si+2));
                    // now Θ refers to x(0)>x(2)
                    fe2_vec[1] = SD::ChengWu::Binarize(fe2_vec[1], x(si+0)-x(si+2))[0];
                    
                    ex xsum = x(si+3)+x(si+4)+x(si+5);
                    //xsum = x(si+5);
                    SD::ChengWu::Projectivize(fe0_vec[0],xs,xsum);
                    SD::ChengWu::Projectivize(fe0_vec[1],xs,xsum);
                    SD::ChengWu::Projectivize(fe2_vec[0],xs,xsum);
                    SD::ChengWu::Projectivize(fe2_vec[1],xs,xsum);
                    
                    return XIntegral(fe0_vec[0]) + XIntegral(fe0_vec[1]) + XIntegral(fe2_vec[0]) + XIntegral(fe2_vec[1]);
                } else return e.map(self);
            })(ret);
            
            return ret;
        } else throw Error("DoPS: only 2-, 3-, 4- massless Phase Space supported.");
    }
    
}

