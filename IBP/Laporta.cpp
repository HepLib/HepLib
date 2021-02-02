/**
 * @file
 * @brief IBP with Leporta
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib::IBP {

    /**
     * @brief export integral to KIRA form
     */
    string Laporta::Fout(const ex & expr) {
        if(!using_uw) {
            ex f = expr;
            if(is_a<lst>(f)) f = F(f);
            else if(expr.match(F(w1,w2))) f = F(f.op(1));
            string fstr = ex2str(f);
            string_replace_all(fstr,"{","");
            string_replace_all(fstr,"}","");
            string_replace_all(fstr,"(","[");
            string_replace_all(fstr,")","]");
            return fstr;
        } else {
            ex idx;
            if(expr.match(F(w1,w2))) idx = expr.op(1);
            else if(expr.match(F(w))) idx = expr.op(0);
            else idx = expr;
            return to_string(i2w[idx]);
        }
    }
    
    /**
     * @brief import integral from KIRA
     */
    ex Laporta::Fin(const string & expr) {
        if(!using_uw) {
            string fstr = expr;
            string_replace_all(fstr,"[","("+to_string(ProblemNumber)+",{");
            string_replace_all(fstr,"]","})");
            return str2ex(fstr);
        } else {
            auto cpos = expr.find("*");
            if(cpos==string::npos) {
                if(expr=="0") return 0;
                throw Error("KIRA::Fin with 0 or * NOT Found.");
            }
            auto wstr = expr.substr(0,cpos);
            unsigned long long weight = stoull(wstr,NULL,0);
            auto oex = str2ex(expr.substr(cpos+1,string::npos));
            return F(ProblemNumber, w2i[weight]) * oex;
        }
    }
    
    /**
     * @brief Export input data for KIRA
     */
    void Laporta::Export() {

        if(Integrals.nops()<1) return;
        
        int pdim = Propagators.nops();
        lst InExternal;
        for(auto ii : Internal) InExternal.append(ii);
        for(auto ii : External) InExternal.append(ii);
        
        if(ISP.nops()<1) {
            for(auto it : Internal) {
                for(auto ii : InExternal) ISP.append(it*ii);
            }
            ISP.sort();
            ISP.unique();
        }
        
        if(ISP.nops() > pdim) {
            cout << "ISP = " << ISP << endl;
            cout << "Propagators = " << Propagators << endl;
            throw Error("UKIRA::Export: #(ISP) > #(Propagators).");
        }
        
        lst sp2s, s2sp, ss;
        int _pic=0;
        for(auto item : ISP) {
            _pic++;
            Symbol si("P"+to_string(_pic));
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst leqns;
        for(int i=0; i<ISP.nops(); i++) { // note NOT pdim
            auto eq = Propagators.op(i).expand().subs(iEpsilon==0); // drop iEpsilon
            eq = eq.subs(sp2s, subs_options::algebraic);
            eq = eq.subs(Replacements, subs_options::algebraic);
            if(eq.has(iWF(w))) throw Error("UKIRA::Export, iWF used in eq.");
            leqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(leqns, ss);
        if(s2p.nops() != ISP.nops()) throw Error("KIRA::Export: lsolve failed.");
       
        if(DSP.nops()<1) {
            for(auto p1 : Internal)
            for(auto p2 : InExternal)
            DSP.append(lst{p1,p2});
        }

        ibps.remove_all(); // no need
        lst nsa;
        for(int i=0; i<pdim; i++) nsa.append(a(i));
        for(auto sp : DSP) {
            auto ilp = sp.op(0);
            auto iep = sp.op(1);
            
            ex ibp = 0;
            symbol ss;
            for(int i=0; i<pdim; i++) {
                auto ns = nsa;
                ns.let_op(i) = nsa.op(i) + 1;
                auto dp = Propagators.op(i).subs(ilp==ss).diff(ss).subs(ss==ilp);
                ibp -= (a(i)+Shift[i]) * F(ns) * dp;
            }
            
            ibp = ibp * iep;
            ibp = ibp.expand();
            ibp = ibp.subs(sp2s, subs_options::algebraic);
            ibp = ibp.subs(Replacements, subs_options::algebraic);
            ibp = ibp.subs(s2p, subs_options::algebraic);
            
            ex res = 0;
            for(int i=0; i<pdim; i++) {
                auto ci = ibp.coeff(iWF(i), 1);
                ci = MapFunction([i](const ex &e, MapFunction &self)->ex {
                    if(e.match(F(w))) {
                        lst tmp = ex_to<lst>(e.op(0));
                        tmp.let_op(i) = tmp.op(i)-1;
                        return F(tmp);
                    } else if(!e.has(F(w))) return e;
                    else return e.map(self);
                })(ci);
                res += ci;
            }
            res += ibp.subs(lst{iWF(w)==0});
            auto cilp = iep.coeff(ilp);
            if(!is_zero(cilp)) res += d*cilp*F(nsa);
            ibps.append(res);
        }
        
cout << ibps << endl;
        

        
    }
    
    /**
     * @brief invoke kira program for reduction
     */
    void Laporta::Run() {
        
    }

    /**
     * @brief import kira result
     */
    void Laporta::Import() {
        
    }

}

