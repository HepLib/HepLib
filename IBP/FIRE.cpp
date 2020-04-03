#include "IBP.h"

namespace HepLib::IBP {

    REGISTER_FUNCTION(a, do_not_evalf_params())
    REGISTER_FUNCTION(F, do_not_evalf_params())

    void FIRE::Reduce() { 
        if(Integrals.nops()<1) return;
        Dimension = Propagators.nops();
        int pn = ProblemNumber;
        int dim = Dimension;
        lst InExternal;
        for(auto ii : Internal) InExternal.append(ii);
        for(auto ii : External) InExternal.append(ii);

        lst sps;
        for(auto it : Internal) {
            for(auto ii : InExternal) sps.append(it*ii);
        }
        sps.sort();
        sps.unique();
        if(sps.nops() != dim) throw Error("FIRE::Reduce: sps failed.");
        
        lst sp2s, s2sp, ss;
        for(auto item : sps) {
            symbol si;
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst eqns, iWFs;
        for(int i=0; i<dim; i++) {
            auto eq = Propagators.op(i).expand();
            eq = eq.subs(sp2s, subs_options::algebraic);
            eq = eq.subs(Replacements, subs_options::algebraic);
            if(eq.has(iWF(w))) throw Error("iWF used in eq.");
            eqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(eqns, ss);
        if(s2p.nops() != dim) throw Error("FIRE::Reduce: lsove failed.");

        IBPs.clear();
        exvector IBPvec;
        lst ns0;
        for(int i=0; i<dim; i++) ns0.append(0);
        for(auto loop : Internal) {
            lst dp_lst;
            for(int i=0; i<dim; i++) {  
                auto s = ex_to<Symbol>(loop);
                dp_lst.append(Propagators.op(i).diff(s));
            }
            
            for(auto iep : InExternal) {
                exmap nc_map;
                for(int i=0; i<dim; i++) { // diff on each propagator
                    auto ns = ns0;
                    ns.let_op(i) = ns.op(i)+1; // note the covention
                    auto tmp = dp_lst.op(i) * iep;
                    tmp = mma_collect(tmp, InExternal);
                    tmp = tmp.subs(Replacements, subs_options::algebraic);
                    tmp = tmp.subs(sp2s, subs_options::algebraic);
                    tmp = tmp.subs(s2p, subs_options::algebraic);
                    tmp = ex(0) - a(i+1)*tmp;

                    for(int j=0; j<dim; j++) {
                        auto cj = tmp.coeff(iWF(j));
                        if(is_zero(cj)) continue;
                        auto cns = ns;
                        cns.let_op(j) = cns.op(j)-1; // note the covention
                        nc_map[cns] = nc_map[cns] + cj;
                    }
                    tmp = tmp.subs(iWF(w)==0); // constant term
                    if(!is_zero(tmp)) nc_map[ns] = nc_map[ns] + tmp;
                }
                
                if(is_zero(loop-iep)) nc_map[ns0] = nc_map[ns0] + d;
                bool ok = false;
                for(auto nc : nc_map) {
                    if(!is_zero(nc.second)) {
                        ok = true;
                        IBPvec.push_back(nc.second);
                    }
                }
                if(ok) IBPs.push_back(nc_map);
            }
        }
        
        auto Variables = gather_symbols(IBPvec);
        
        ostringstream start;
        start << "Null" << endl << endl;
        start << "ExampleDimension[" << pn << "]=" << Internal.nops() * InExternal.nops() << endl << endl;
        start << "ProblemNumber=" << pn << endl << endl;
        
        // .start - SBasisL
        PermutationsR(2, dim, [dim,pn,&start](const int *ns) {
            start << "SBasisL[" << pn << ",{";
            for(int i=0; i<dim; i++) start << (ns[i]<1 ? -1 : 1) << (i<dim-1 ? "," : "");
            start << "}]=0" << endl << endl;
        });
        
        // .start - SBasis0L, SBasis0D && SBasis0C
        start << "SBasis0L[" << pn << "]=" << IBPs.size() << endl << endl;
        ostringstream oss;
        for(int i=0; i<IBPs.size(); i++) {
            start << "SBasis0D[" << pn << "," << (i+1) << "]=";
            lst items;
            for(auto kv : IBPs[i]) {
                items.append(kv.first);
                oss << "SBasis0C[" << pn << "," << (i+1) << "," << kv.first << "]=" << kv.second << endl << endl;
            }
            start << items << endl << endl;
        }
        start << oss.str();
        
        // .start - SBasisS
        start << "SBasisS[" << pn << "]={{{";
        // 1,2,3,...
        for(int i=0; i<dim; i++) start << (i+1) << (i<dim-1 ? "," : "");
        start << "},{";
        // 1,1,1,...
        for(int i=0; i<dim; i++) start << 1 << (i<dim-1 ? "," : "");
        start << "},{";
        // 0,0,0,...
        for(int i=0; i<dim; i++) start << 0 << (i<dim-1 ? "," : "");
        start << "}}}" << endl << endl;
        
        // .start - SBasisR
        start << "SBasisR[" << pn << ",{";
        for(int i=0; i<dim; i++) start << "-1" << (i<dim-1 ? "," : "");
        start << "}]=True" << endl << endl;
        
        // .start - Others
        start << "SBasisRL[" << pn << "]=0" << endl << endl;
        start << "HPI[" << pn << "]={}" << endl << endl;
        
        string sss = start.str();
        string_replace_all(sss, "=", " = ");
        string_replace_all(sss, ",", ", ");
        
        string spn = to_string(pn);
        
        // .config
        ostringstream config;
        config << "#threads 4" << endl;
        config << "#fermat fer64" << endl;
        config << "#variables ";
        bool first = true;
        for(auto v : Variables) { config << (first ? "" : ",") << v; first=false; }
        config << endl;
        config << "#database db" << pn << endl;
        config << "#bucket 20" << endl;
        config << "#start" << endl;
        config << "#problem " << pn << " " << pn << ".start" << endl;
        config << "#integrals " << pn << ".intg" << endl;
        config << "#output " << pn << ".tables" << endl;
        
        // *.intg
        ostringstream intg;
        intg << "{";
        for(int i=0; i<Integrals.nops(); i++) {
            intg << "{" << pn << "," << Integrals[i] << (i<Integrals.nops()-1 ? "}," : "}");
        }
        intg << "}" << endl;
        
        if(WorkingDir.length()<1) WorkingDir = to_string(getpid());
        system(("mkdir -p "+WorkingDir).c_str());
        ofstream start_out(WorkingDir+"/"+spn+".start");
        start_out << sss << endl;
        start_out.close();
        
        ofstream config_out(WorkingDir+"/"+spn+".config");
        config_out << config.str() << endl;
        config_out.close();
        
        ofstream intg_out(WorkingDir+"/"+spn+".intg");
        intg_out << intg.str() << endl;
        intg_out.close();
        
        ostringstream cmd;
        cmd << "cd " << WorkingDir << " && FIRE5 -c " << pn << " >/dev/null";
        system(cmd.str().c_str());
        system(("rm -rf "+WorkingDir+"/db"+to_string(pn)).c_str());
        
        ifstream ifs(WorkingDir+"/"+spn+".tables");
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        
        string_replace_all(ostr, "\"", "");
        
        Parser tp;
        auto tp_lst = tp.Read(ostr);
        exmap id2F;
        for(auto item : tp_lst.op(1)) {
            id2F[item.op(0)] = F(item.op(1).op(0), item.op(1).op(1));
        }
        
        for(auto item : tp_lst.op(0)) {
            ex left = item.op(0).subs(id2F);
            ex right = 0;
            for(auto it : item.op(1)) {
                right += it.op(0).subs(id2F) * it.op(1);
            }
            if(is_zero(left-right)) MasterIntegrals.append(left);
            else Rules.append(left==right);
        }
    }
    
    ex FIRE::UF(ex idx) {
        ex ft = 0;
        int nps = Propagators.nops();
        for(int i=0; i<nps; i++) {
            if(is_zero(idx.op(i))) continue;
            ft -= (is_zero(idx.op(i)-1) ? ex(1) : iWF(idx.op(i))) * x(i) * Propagators.op(i);
        }

        ex ut = 1;
        for(int i=0; i<Internal.nops(); i++) {
            ft = ft.expand();
            auto t2 = ft.coeff(Internal.op(i), 2);
            auto t1 = ft.coeff(Internal.op(i), 1);
            auto t0 = ft.subs(Internal.op(i)==0);
            ut *= t2;
            if(is_zero(t2)) return lst{0,0};
            ft = normal(t0-t1*t1/(4*t2));
        }
        ft = subs(ex(0)-ut*ft, Replacements, subs_options::algebraic);
        ft = normal(ft);
        ut = subs(ut, Replacements, subs_options::algebraic);
        ut = normal(ut);
        ex uf = ut*ft;

        ex_is_less comp;
        Permutations(nps, [comp,nps,&uf,&ft,&ut](const int *ns) {
            exmap x2x;
            for(int i=0; i<nps; i++) x2x[x(i)]=x(ns[i]);
            ex uf2 = uf.subs(x2x);
            if(comp(uf2, uf)) {
                uf = uf2;
                ft = ft.subs(x2x);
                ut = ut.subs(x2x);
            }
        });

        return lst{ut, ft};
    }  
    
    lst FIRE::FindRules(const vector<FIRE> & fs, bool mi) {
        map<ex,lst,ex_is_less> group;;
        if(mi) {
            for(auto fi : fs) {
                for(auto mi : fi.MasterIntegrals) group[fi.UF(mi.subs(F(w1,w2)==w2))].append(mi);
            }
        } else {
            for(auto fi : fs) {
                for(auto idx : fi.Integrals) group[fi.UF(idx)].append(F(fi.ProblemNumber,idx));
            }
        }
        
        lst rules;
        for(auto g : group) {
            lst gs = ex_to<lst>(g.second);
            if(gs.nops()<2) continue;
            for(int i=1; i<gs.nops(); i++) rules.append(gs.op(i)==gs.op(0));
        }
        return rules;
    }
    
    lst FIRE::FindRules(const vector<FIRE*> & fs, bool mi) {
        map<ex,lst,ex_is_less> group;;
        if(mi) {
            for(auto fi : fs) {
                for(auto mi : fi->MasterIntegrals) group[fi->UF(mi.subs(F(w1,w2)==w2))].append(mi);
            }
        } else {
            for(auto fi : fs) {
                for(auto idx : fi->Integrals) group[fi->UF(idx)].append(F(fi->ProblemNumber,idx));
            }
        }
        
        lst rules;
        for(auto g : group) {
            lst gs = ex_to<lst>(g.second);
            if(gs.nops()<2) continue;
            for(int i=1; i<gs.nops(); i++) rules.append(gs.op(i)==gs.op(0));
        }
        return rules;
    }

}
