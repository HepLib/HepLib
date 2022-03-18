/**
 * @file
 * @brief Basic Functions, extend GiNaC
 */

#include "DE.h"

namespace HepLib {
        
    DESS::DESS(IBP & _ibp, const symbol & _x) : ibp(_ibp), x(_x), DE(_x) { } 
    
    void DESS::InitDE() {
        if(Verbose>1) cout << "  \\--Generating DE @ " << now() << endl;
        
        if(ibp.MIntegrals.nops()<1) {
            ibp.Reduce();
            ibp.RM(true); // keep .start & .config
        }
        ibp.FindRules(true);
        Rules = ibp.Rules;
        lst ns;
        for(int i=0; i<ibp.Propagators.nops(); i++) ns.append(a(i));
        auto ibpD = ibp.D(x, ns);
        while(true) {
            ibp.Rules.remove_all();
            MIntegrals = ibp.MIntegrals;
            sort_lst(MIntegrals);
            ibp.MIntegrals.remove_all();
            ibp.Integrals.remove_all();
            lst dmis;
            for(auto mi : MIntegrals) {
                auto nsi = ex_to<lst>(mi.op(1));
                ex dmi = ibpD.subs(ns,nsi,nopat);
                dmis.append(dmi);
            }
            exset fs;
            find(dmis, F(w1,w2), fs);
            for(auto fi : fs) ibp.Integrals.append(fi.op(1));
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
        if(Verbose>1) cout << "  \\--Generated DE @ " << now() << endl;
    }
        
}

