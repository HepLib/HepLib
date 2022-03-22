/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {
        
    AMF::AMF(IBP & _ibp) : ibp(_ibp), x(iet) { } 
    
    void AMF::InitDE() {
        if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--Generating DE @ " << now() << endl;
        
        if(ibp.MIntegrals.nops()<1) {
            ibp.Reduce();
            ibp.RM(false); // since Propagators wiil be changed
        }
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
            if(!In_GiNaC_Parallel && Verbose>1) cout << "  \\--DE Poles: " << roots << endl;
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
        DE de(x, Mat);
        de.d0 = d0;
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
        oo.d0 = d0;
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
        auto ooMM = oo.Series(xoo, xN);
        auto ooT = oo.MatT();
        auto ooU0 = ooMM.second;
        
        // ooC
        matrix ooC; // F in ooC should be treated as Bubble type
        if(true) {
            matrix MIs(matN,1);
            for(int i=0; i<matN; i++) MIs(i,0) = dlst.op(i) * MIntegrals.op(i);
            if(is_a<numeric>(d0)) ooC = ooT.mul(ooU0).inverse().mul(MIs);
            else ooC = fermat_inv(ooT.mul(ooU0)).mul(MIs);
            xpow(ooC,x);
            HepLib::subs(ooC,x==0,nopat);
        }
    
        // J(xoo)
        matrix ooTUC = ooMM.first.mul(ooC);
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
        o.d0 = d0;
        o.WDigits = WDigits;
        ex xo = pts.op(0);
        auto oMM = o.Series(xo, xN);
        matrix oT = o.MatT();
        matrix oU0 = oMM.second;
        matrix ioT = oT.inverse();
        if(WDigits>0) ioT = ex_to<matrix>(subs(ioT,x==xo.evalf()));
        else ioT = ex_to<matrix>(subs(ioT,x==xo));
        matrix ioTU = oMM.first.inverse().mul(ioT);

        // BC at origin
        matrix oC = ioTU.mul(MatU.mul(ooTUC));
        
        // only pick up x^integer parts
        for(int i=0; i<oU0.nops(); i++) {
            ex item = oU0.op(i);
            if(!is_a<mul>(item)) item = lst{item};
            for(auto it : item) {
                if(it.match(pow(x,w1)) && !it.op(1).info(info_flags::integer)) {
                    oU0.let_op(i) = 0;
                    goto next_U0;
                }
            }
            next_U0: ;
        }
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
    
    lst AMF::Evaluate(const lst & d0s, bool parallel) {
        lst res;
        if(parallel && d0s.nops()>1) {
            auto od0 = d0;
            auto res_vec = GiNaC_Parallel(d0s.nops(), [&](int idx)->ex {
                d0 = d0s.op(idx);
                return Evaluate();
            }, "AMF");
            d0 = od0;
            res = vec2lst(res_vec);
        } else {
            auto od0 = d0;
            for(auto di : d0s) {
                d0 = di;
                res.append(Evaluate());
            }
            d0 = od0;
        }
        return res;
    }
        
}

