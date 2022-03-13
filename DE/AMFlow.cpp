/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {
        
    AMFlow::AMFlow(IBP & _ibp) : ibp(_ibp), o(iet), oo(iet) { } 
    
    void AMFlow::InitDE() {
        if(ibp.MIntegrals.nops()<1) ibp.Reduce();
        ibp.FindRules(true);
        Rules = ibp.Rules;
        int matN;
        for(int i=0; i<ibp.Propagators.nops(); i++) {
            ibp.Propagators.let_op(i) = ibp.Propagators.op(i) + iet;
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
            system(("rm -rf "+ibp.WorkingDir).c_str());
            ibp.Reduce();
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
                o.Mat = oo.Mat = mat;
                system(("rm -rf "+ibp.WorkingDir).c_str());
                break;
            }
        }
    }
    
    matrix AMFlow::FSS(const int xn) {
        if(o.Mat.nops()<1) InitDE();
        o.Fuchsify();
        o.Shear();
        return o.Series(o.x, xn);
    }
    
    void AMFlow::Scale() {
        lst diag;
        int nloop = ibp.Internal.nops();
        for(auto mi : MIntegrals) {
            ex v = 0;
            for(auto ni : mi.op(1)) v += ni;
            diag.append(pow(oo.x,nloop*d/2-v));
        }
        //oo.Apply(diag);
        oo.x2y(1/oo.x);
        oo.xpow();
        oo.Fuchsify();
        oo.Shear();
        oo.info();
        auto m = oo.Series(o.x,1);
        cout << m << endl << endl;
        cout << oo.MatT() << endl << endl;
        cout << MIntegrals << endl;
    }
        
}

