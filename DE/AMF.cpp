/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {
        
    AMF::AMF(IBP & _ibp) : ibp(_ibp), x(iet) { } 
    
    void AMF::InitDE() {
        if(Verbose>1) cout << "  \\--Generating DE @ " << now() << endl;
        
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
                Mat = ex_to<matrix>(subs(mat,d==d0,nopat));
                break;
            }
        }
        system(("rm -rf "+ibp.WorkingDir).c_str());
        if(Verbose>1) cout << "  \\--Generated DE @ " << now() << endl;
        
        if(Verbose>10) cout << "  \\--DE Poles starting @ " << now() << endl;
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
                if (deg == 0) { }
                else if (deg == 1) {
                    ex c0 = b.coeff(x, 0);
                    ex c1 = b.coeff(x, 1);
                    roots.insert(normal(-c0/c1));
                } else throw Error("AMF::InitDE, higher powers found.");
            }
            if(Verbose>1) cout << "  \\--DE Poles: " << roots << endl;
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
        if(Verbose>10) {
            cout << "  \\--AMF Points: " << NN(pts,2) << endl;
            cout << "  \\--DE Poles finished @ " << now() << endl;
        }
    }
    
    matrix AMF::RU(const ex & x1, const ex & x2) {
        ex dis = x1-x2;
        DE de(x, Mat);
        de.Precision = Precision;
        de.x2y(x+x2); // transform x2 to 0
        auto mat = de.Taylor(dis,xN);
        return mat;
    }
    
    lst AMF::Evaluate() {
        int matN = Mat.rows();
        DE oo(Mat,x);
        oo.x2y(1/x); // infinity to origin
        oo.Precision = Precision;
        oo.d0 = d0;
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
        
        // make sure U0 is diagonal
        for(int r=0; r<matN; r++) for(int c=0; c<matN; c++) if(r!=c && !is_zero(ooU0(r,c))) 
            throw Error("AMF::oo2o, U0 at infinity is NOT diagonal");
            
        auto Ti = ooT.inverse();
        matrix ooC(matN, 1); 
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
            ci = oo.xpow(ci/ooU0(i,i)).subs(x==0, nopat);
            ooC(i,0) = ci; // F in ooC should be treated as Bubble type
        }
        
        // J(xoo)
        matrix ooTUC = ex_to<matrix>(subs(ooT, x==xoo)).mul(ooMM.first).mul(ooC);
        
        // Middle U matrix
        matrix MatU = ex_to<matrix>(unit_matrix(matN));
        for(int i=0; i<npts-1; i++) MatU = MatU.mul(RU(pts.op(i),pts.op(i+1)));
        
        DE o(Mat,x);
        o.Precision = Precision;
        ex xo = pts.op(0);
        auto oMM = o.Series(xo, xN);
        auto oT = o.MatT();
        auto oU0 = oMM.second;
        auto oTU = ex_to<matrix>(subs(oT, x==xo)).mul(oMM.first);
        if(Precision>0) oTU = oTU.inverse();
        else oTU = fermat_inv(oTU);

        // BC at origin
        matrix oC = oTU.mul(MatU).mul(ooTUC);
        
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
        
        lst res; // result for master integrals
        oC = ex_to<matrix>(oT.mul(oU0)).mul(oC);
        for(int i=0; i<matN; i++) res.append(oC(i,0));
        return res;
    }
        
}

