/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"
#include <cmath>

namespace HepLib {
    
    using namespace EoD;
        
    AMF::AMF(IBP & _ibp) : ibp(_ibp), x(iet) { } 
    
    void AMF::InitDE() {
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Generating DE @ " << now() << endl;
        
        if(ibp.MIntegrals.nops()<1) {
            ibp.Reduce();
            ibp.RM(false); // since Propagators will be changed
        }
        if(ibp.MIntegrals.nops()<1) throw Error("No MI Found, maybe No need DE!");
        ibp.FindRules(true);
        Rules = ibp.Rules;
        for(int i=0; i<ibp.Propagators.nops(); i++) {
            ibp.Propagators.let_op(i) = ibp.Propagators.op(i) + x;
        }
        while(true) {
            ibp.Rules.remove_all();
            MIntegrals = ibp.MIntegrals;
            sort_lst(MIntegrals);
            ibp.MIntegrals.remove_all();
            ibp.Integrals.remove_all();
            lst dmis;
            for(auto mi : MIntegrals) {
                ex dmi = 0;
                auto ns0 = mi.op(1);
                for(int i=0; i<ns0.nops(); i++) {
                    if(is_zero(ns0.op(i))) continue;
                    auto ns = ns0;
                    ns.let_op(i) = ns0.op(i)+1;
                    dmi -= ns0.op(i)*F(mi.op(0),ns);
                    ibp.Integrals.append(ns);
                }
                dmis.append(dmi);
            }
            ibp.Reduce();
            ibp.RM(true); // keep .start & .config
            ibp.FindRules(true);
            sort_lst(ibp.MIntegrals);

            if(ibp.MIntegrals==MIntegrals) {
                int matN = MIntegrals.nops();
                matrix mat(matN,matN);
                for(int r=0; r<matN; r++) {
                    auto dmi = dmis.op(r).subs(ibp.Rules,nopat);
                    dmi = collect_ex(dmi, F(w1,w2));
                    ex chk = 0;
                    for(int c=0; c<matN; c++) {
                        mat(r,c) = dmi.coeff(MIntegrals.op(c));
                        chk += mat(r,c) * MIntegrals.op(c);
                    }
                    if(!is_zero(normal(chk-dmi))) throw Error("AMFlow::InitDE, Check failed");
                }
                Mat = mat;
                matrix_map_inplace(Mat, [](const ex & e) { return normal(e); });
                break;
            }
        }
        system(("rm -rf "+ibp.WorkingDir).c_str());
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Generated DE @ " << now() << endl;
        
        if(!In_GiNaC_Parallel && Verbose>10) cout << "  \\--DE Poles starting @ " << now() << endl;
        if(true) { // get poles of Mat
            ex den = matrix_den_lcm(Mat);
            exvector fvec;
            if (is_a<mul>(den)) for (const auto &f : den) fvec.push_back(f);
            else fvec.push_back(den);
            exset roots;
            for (const auto &f : fvec) {
                ex b = f; 
                int n = 1;
                if (is_a<power>(f)) {
                    b = f.op(0);
                    n = ex_to<numeric>(f.op(1)).to_int();
                } 
                int deg = b.degree(x);
                if (deg==0) { }
                else if (deg==1) {
                    ex c0 = b.coeff(x, 0);
                    ex c1 = b.coeff(x, 1);
                    roots.insert(normal(-c0/c1));
                } else if(deg==2) {
                    ex c0 = b.coeff(x, 0);
                    ex c1 = b.coeff(x, 1);
                    ex c2 = b.coeff(x, 2);
                    roots.insert(Rationalize(evalf((-c1+sqrt(c1*c1-4*c2*c0))/(2*c2)),10));
                    roots.insert(Rationalize(evalf((-c1-sqrt(c1*c1-4*c2*c0))/(2*c2)),10));
                } else {
                    cout << "current factor: " << f << endl;
                    throw Error("AMF::InitDE, higher powers found.");
                }
            }
            lst rs;
            for(auto ri : roots) rs.append(ri);
            sort_lst(rs);
            if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--DE Poles: " << NN(rs,2) << endl;
            pts.remove_all();
            if(roots.size()>0) {
                ex max = -1, min = -1;
                for(auto r : roots) {
                    if(is_zero(r)) continue;
                    ex ar = abs(r);
                    if(max<0 || ar>max) max = ar;
                    if(min<0 || ar<min) min = ar;
                }
                pts.remove_all();
                ex x0 = I*2*max; // the last point
                x0 = Rationalize(x0, 20);
                pts.append(x0);
                
                while(true) {
                    ex mm = -1;
                    for(auto r : roots) {
                        ex ar = abs(r-x0);
                        if(mm<0 || ar<mm) mm = ar;
                    }
                    x0 -= I*mm/2;
                    x0 = Rationalize(x0, 20);
                    pts.prepend(x0);
                    if(x0/I<0) throw Error("AMF::InitDE, et<0 FOUND.");
                    if(abs(x0)<min) break;
                }
            } else throw Error("AMF::InitDE, NO root found.");
        }
        if(!In_GiNaC_Parallel && Verbose>10) {
            cout << "  \\--AMF Points: " << NN(pts,2) << endl;
            cout << "  \\--DE Poles finished @ " << now() << endl;
        }
    }
    
    matrix AMF::RU(const ex & x1, const ex & x2) {
        ex dis = x1-x2;
        DE de(Mat,x);
        if(d!=d0) de.subs(d==d0, nopat);
        de.WDigits = WDigits;
        auto mat = de.Taylor(x2,dis,xN);
        return mat;
    }
    
    lst AMF::Evaluate() {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        int matN = Mat.rows();
        
        // DE at infinity
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ infinity ..." << RESET << endl;
        DE oo(Mat,x);
        if(d!=d0) oo.subs(d==d0, nopat);
        oo.WDigits = WDigits;
        oo.x2y(1/x); // infinity to origin
        lst dlst;
        if(true) { // rescale MI
            int dim = ibp.Propagators.nops();
            int nloop = ibp.Internal.nops();
            for(int i=0; i<matN; i++) {
                ex v = -nloop * d0 / 2;
                ex ns = MIntegrals.op(i).op(1);
                for(int j=0; j<dim; j++) v += ns.op(j);
                dlst.append(pow(x,v));
            }
            oo.Apply(dlst);
        }
        oo.xpow();
        int npts = pts.nops();
        ex xoo = 1/pts.op(npts-1);
        auto ooU = oo.Series(xoo, xN);
        auto ooT = oo.MatT();
        auto ooU0 = oo.Series(x,0);       

        // ooC 
        matrix ooC; // F in ooC should be treated as Bubble type
        int ooCver = 2;
        // make sure U0 is diagonal
        if(false)
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) if(r!=c && !is_zero(ooU0(r,c))) 
            throw Error("AMF::oo2o, U0 at infinity is NOT diagonal");
            
        if(ooCver==1) {
            auto Ti = ooT.inverse();
            ooC = matrix(matN, 1); 
            for(int i=0; i<matN; i++) {
                if(true) { // check ooU0
                    if(!ooU0(i,i).is_equal(1)) {
                        if(!ooU0(i,i).match(pow(x,w))) {
                            cout << endl << "ooU0 = " << ooU0 << endl;
                            throw Error("AMF::oo2o, wrong pattern found.");
                        }
                        if(ooU0(i,i).op(1)>0) {
                            cout << endl << "ooU0 = " << ooU0 << endl;
                            throw Error("AMF::oo2o, wrong pattern found.");
                        }
                    }
                }
                ex ci = 0;
                for(int j=0; j<matN; j++) ci += Ti(i,j) * dlst.op(j) * MIntegrals.op(j);
                ci = xpow(ci/ooU0(i,i),x).subs(x==0, nopat);
                ooC(i,0) = ci; // F in ooC should be treated as Bubble type
            }
        } else if(ooCver==2) {
            matrix MIs(matN,1);
            for(int i=0; i<matN; i++) MIs(i,0) = dlst.op(i) * MIntegrals.op(i);
            if(is_a<numeric>(d0)) ooC = ooT.mul(ooU0).inverse().mul(MIs);
            else ooC = fermat_inv(ooT.mul(ooU0)).mul(MIs);
            xpow(ooC,x);
            HepLib::subs(ooC,x==0,nopat);
        }
    
        // J(xoo)
        matrix ooTUC = ooU.mul(ooC);
        if(WDigits>0) ooTUC = ex_to<matrix>(subs(ooT, x==xoo.evalf())).mul(ooTUC);
        else ooTUC = ex_to<matrix>(subs(ooT, x==xoo)).mul(ooTUC);
        
        // Middle U matrix
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ regular points ..." << RESET << endl;
        matrix MatU = ex_to<matrix>(unit_matrix(matN));
        for(int i=0; i<npts-1; i++) MatU = MatU.mul(RU(pts.op(i),pts.op(i+1)));
                
        // DE at origin
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "AMF @ origin ..." << RESET << endl;
        DE o(Mat,x);
        if(d!=d0) o.subs(d==d0, nopat);
        o.WDigits = WDigits;
        ex xo = pts.op(0);
        auto oU = o.Series(xo, xN);
        matrix oT = o.MatT();
        matrix ioT = oT.inverse();
        if(WDigits>0) ioT = ex_to<matrix>(subs(ioT,x==xo.evalf()));
        else ioT = ex_to<matrix>(subs(ioT,x==xo));
        matrix ioTU = oU.inverse().mul(ioT);
        
        // Final C at origin
        matrix oC = ioTU.mul(MatU.mul(ooTUC));
        
        // take x->0 limit at origin
        auto pr = prank(oT,x);
        pr++;
        if(pr<0) pr = 0;
        for(auto ev : o.EigenValues()) {
            if(!ev.info(info_flags::integer)) continue;
            pr -= ex_to<numeric>(ev).to_int();
            break;
        }
        matrix oU0 = o.Series(x,pr,lst{0}); // only pick up x^integer
        auto oTU = oT.mul(oU0);
        xpow(oTU,x);
        if(WDigits>0) oTU = ex_to<matrix>(oTU.evalf());
        HepLib::subs(oTU,x==0,nopat);
        
        lst res; // result for master integrals
        oC = oTU.mul(oC);
        for(int i=0; i<matN; i++) res.append(oC(i,0));
        Digits = oDigits;
        return res;
    }
    
    lst AMF::FitEps(const lst & eps, int lp, bool parallel) {
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        exvector eps_vec(eps.begin(), eps.end());
        int nmi = MIntegrals.nops();
        exvector mis_vec[nmi];
        if(parallel && eps.nops()>1) {
            auto od0 = d0;
            auto res_vec = GiNaC_Parallel(eps.nops(), [&](int idx)->ex {
                d0 = 4-2*eps.op(idx);
                return Evaluate();
            }, "AMF");
            d0 = od0;
            for(auto mis : res_vec) {
                for(int i=0; i<nmi; i++) mis_vec[i].push_back(mis.op(i));
            }
        } else {
            auto od0 = d0;
            for(auto epi : eps) {
                d0 = 4-2*epi;
                if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--" << WHITE << "ep = " << epi << RESET << endl;
                auto mis = Evaluate();
                for(int i=0; i<nmi; i++) mis_vec[i].push_back(mis.op(i));
            }
            d0 = od0;
        }
        
        if(!In_GiNaC_Parallel && Verbose>0) cout << "  \\--Final PolynomialFit ..." << endl;
        if(WDigits>0) Digits = WDigits;
        lst mis_lst;
        for(int i=0; i<nmi; i++) {
            int tn = eps.nops()-1;
            auto cs = PolynomialFit(eps_vec, mis_vec[i], tn, lp);
            ex mi = 0;
            for(int i=0; i<tn; i++) mi += pow(ep,lp+i) * cs.op(i);
            mis_lst.append(mi);
        }
        Digits = oDigits;
        return mis_lst;
    }
    
    lst AMF::FitEps(int goal, int order, bool parallel) { // form AMFlow
        auto oDigits = Digits;
        if(WDigits>0) Digits = WDigits;
        int nloop = ibp.Internal.nops();
        numeric nn = 5*order/numeric(2)+2*nloop;
        int n = ceil(nn.to_double());
        if(n>100) throw Error("FitEps: too large order.");
        ex eps0 = pow(10, -nloop/numeric(2)-goal/numeric(order+1));
        eps0 = Rationalize(evalf(eps0));
        int lp = -2*nloop;
        lst eps;
        for(int i=0; i<n; i++) eps.append(eps0 + eps0*ex(i+1)/100);
        numeric nsp = (n+2*nloop)*(nloop/numeric(2)+goal/numeric(order+1));
        int sp = ceil(nsp.to_double());
        if(sp<30) sp = 30;
        WDigits = 2*sp;
        xN = 4*sp;
        auto res = FitEps(eps,lp,parallel);
        Digits = oDigits;
        return res;
    }  
        
}

