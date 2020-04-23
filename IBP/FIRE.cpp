/**
 * @file
 * @brief IBP with FIRE
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-20
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib::IBP {

    static ex f_expand(const ex & arg1, const ex & arg2, unsigned options) {
        return F(arg1,arg2).hold();
    }

    static void a_print(const ex & ex_in, const print_context & c) {
        c.s << "a[" << ex_in << "]";
    }

    REGISTER_FUNCTION(a, do_not_evalf_params().print_func<print_dflt>(a_print))
    REGISTER_FUNCTION(F, do_not_evalf_params().expand_func(f_expand))

    /**
     * @brief Do IBP Reduction 
     */
    void FIRE::Reduce() { 
        Dimension = Propagators.nops();
        if(Integrals.nops()<1) return;
        int pn = 0; // to avoid unsigned short overflow in FIRE
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
        if(sps.nops()<dim) {
            lst sps_ext;
            for(auto it : External) {
                for(auto ii : External) sps_ext.append(it*ii);
            }
            sps_ext.sort();
            sps_ext.unique();
            for(auto item : sps_ext) {
                auto item2 = item.subs(Replacements, subs_options::algebraic);
                if(is_zero(item-item2)) sps.append(item);
            }
            sps.sort();
            sps.unique();
        }
        if(sps.nops() != dim) {
            cout << "sps = " << sps << endl;
            cout << "Propagators = " << Propagators << endl;
            throw Error("FIRE::Reduce: sps failed.");
        }
        
        lst sp2s, s2sp, ss;
        for(auto item : sps) {
            symbol si;
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst eqns, iWFs;
        for(int i=0; i<dim; i++) {
            auto eq = Propagators.op(i).expand().subs(iEpsilon==0); // drop iEpsilon
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
                    tmp = tmp.subs(D==d); // replace D to d

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

        start << "ExampleDimension[" << pn << "]=" << dim << endl << endl;
        start << "ProblemNumber=" << pn << endl << endl;
        
        // .start - SBasisL
        if(Version==5) {
            PermutationsR(2, dim, [dim,pn,&start](const int *ns) {
                start << "SBasisL[" << pn << ",{";
                for(int i=0; i<dim; i++) start << (ns[i]<1 ? -1 : 1) << (i<dim-1 ? "," : "");
                start << "}]=0" << endl << endl;
            });
        }
        
        // .start - SBasis0L, SBasis0D && SBasis0C
        start << "SBasis0L[" << pn << "]=" << IBPs.size() << endl << endl;
        ostringstream oss;
        for(int i=0; i<IBPs.size(); i++) {
            start << "SBasis0D[" << pn << "," << (i+1) << "]=";
            lst items;
            for(auto kv : IBPs[i]) {
                items.append(kv.first);
                if(Version==5) {
                    oss << "SBasis0C[" << pn << "," << (i+1) << "," << kv.first << "]=" << 
                    collect_common_factors(kv.second.normal()) << endl << endl;
                } else {
                    lst olst;
                    auto tmp = mma_collect(kv.second, a(w), true, true);
                    if(!is_a<add>(tmp)) tmp = lst{tmp};
                    for(auto item : tmp) {
                        auto cc = item.subs(lst{coVF(w)==1, coCF(w)==w});
                        auto cv = item.subs(lst{coVF(w)==w, coCF(w)==1});
                        cc = collect_common_factors(cc.normal());
                        if(is_zero(cc)) continue;
                        if(is_zero(cv-1)) cv=0;
                        else cv = cv.subs(a(w)==w);
                        olst.append(lst{cc, cv});
                    }
                    oss << "SBasis0C[" << pn << "," << (i+1) << "," << kv.first << "]=" << olst << endl << endl;
                }
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
        lst Rlst;
        Rlst.append(lst{});
        for(int i=0; i<dim; i++) {
            let_op_append(Rlst, 0, -1);
        }
        for(auto lpi : Internal) {
            vector<int> ns_vec;
            lst ns0;
            for(int i=0; i<dim; i++) ns0.append(1);
            for(int i=0; i<dim; i++) {
                if(Propagators.op(i).has(lpi)) ns0.let_op(i) = -1;
                else ns_vec.push_back(i);
            }
            for(int n=0; n<std::pow(2,ns_vec.size()); n++) {
                int cn = n;
                lst ns1 = ns0;
                for(int j=0; j<ns_vec.size(); j++) {
                    if((cn%2)==1) ns1.let_op(ns_vec[j]) = -1;
                    cn /= 2;
                }
                Rlst.append(ns1);
            } 
        }
        
        // handle Cut Propagators
        if(Cuts.nops()>0) {
            //Rlst.remove_all();  // TODO: check this one
            for(auto cx : Cuts) {
                int ci = ex_to<numeric>(cx-1).to_int(); // start from 1 in Cuts
                lst ns0;
                for(int i=0; i<dim; i++) ns0.append(1);
                ns0.let_op(ci) = -1;
                for(int n=0; n<std::pow(2,dim-1); n++) {
                    int cn = n;
                    lst ns1 = ns0;
                    for(int j=0; j<dim; j++) {
                        if(ci==j) continue;
                        if((cn%2)==1) ns1.let_op(j) = -1;
                        cn /= 2;
                    }
                    Rlst.append(ns1);
                } 
            }
        }
        
        Rlst.sort();
        Rlst.unique();
        
        for(auto iR : Rlst) {
            start << "SBasisR[" << pn << ",{";
            for(int i=0; i<dim; i++) start << iR.op(i) << (i<dim-1 ? "," : "");
            start << "}]=True" << endl << endl;
        }
        
        // .start - Others
        start << "SBasisRL[" << pn << "]=0" << endl << endl;
        start << "HPI[" << pn << "]={}" << endl << endl;
        
        string sss = start.str();
        string_replace_all(sss, "=", " = ");
        string_replace_all(sss, ",", ", ");
        
        // .config
        ostringstream config;
        if(Version>5) config << "#compressor none" << endl;
        config << "#threads 1" << endl;
        config << "#fthreads 4" << endl;
        //config << "#fermat fer64" << endl;
        config << "#variables ";
        bool first = true;
        for(auto v : Variables) { config << (first ? "" : ",") << v; first=false; }
        config << endl;
        config << "#database db" << ProblemNumber << endl;
        if(Version>5) config << "#pos_pref 5" << endl;
        if(Version==5) config << "#bucket 20" << endl;
        config << "#start" << endl;
        config << "#problem " << pn << " " << ProblemNumber << ".start" << endl;
        config << "#integrals " << ProblemNumber << ".intg" << endl;
        config << "#output " << ProblemNumber << ".tables" << endl;
        
        // *.intg
        ostringstream intg;
        intg << "{";
        for(int i=0; i<Integrals.nops(); i++) {
            intg << "{" << pn << "," << Integrals[i] << (i<Integrals.nops()-1 ? "}," : "}");
        }
        intg << "}" << endl;
        
        string spn = to_string(ProblemNumber);
        
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
        cmd << "cd " << WorkingDir << " && $(which FIRE" << Version << ") -c " << ProblemNumber << " >/dev/null";
        system(cmd.str().c_str());
        system(("rm -rf "+WorkingDir+"/db"+to_string(ProblemNumber)).c_str());
        
        ifstream ifs(WorkingDir+"/"+spn+".tables");
        string ostr((istreambuf_iterator<char>(ifs)), (istreambuf_iterator<char>()));
        ifs.close();
        
        string_replace_all(ostr, "\"", "");
        
        Parser tp;
        auto tp_lst = tp.Read(ostr);
        exmap id2F;
        for(auto item : tp_lst.op(1)) {
            if(!is_zero(item.op(1).op(0))) throw Error("FIRE::Reduce: pn is NOT 0.");
            id2F[item.op(0)] = F(ProblemNumber, item.op(1).op(1));
        }
        
        for(auto item : tp_lst.op(0)) {
            ex left = item.op(0).subs(id2F);
            ex right = 0;
            for(auto it : item.op(1)) {
                right += it.op(0).subs(id2F) * it.op(1);
            }
            right = right.subs(d==D); // replace d back to D
            if(is_zero(left-right)) MasterIntegrals.append(left);
            else Rules.append(left==right);
        }
        MasterIntegrals.sort();
        MasterIntegrals.unique();
        
    }
    
    /**
     * @brief UF function, from FIRE.m
     * @param idx exponent for the internal Propagator
     * @return lst of {U, F}
     */
    ex FIRE::UF(const ex & idx) const {
        ex ft = 0;
        int nps = Propagators.nops();
        map<int, vector<int>> ngrp;
        for(int i=0; i<nps; i++) ngrp[ex_to<numeric>(idx.op(i)).to_int()].push_back(i);
        
        int nxi=0;
        map<int, vector<int>> pgrp;
        for(auto kv : ngrp) {
            if(kv.first==0) continue;
            for(auto i : kv.second) {
                if(!is_zero(idx.op(i)-1)) ft -= iWF(idx.op(i)) * x(nxi) * Propagators.op(i);
                else ft -= x(nxi) * Propagators.op(i);
                pgrp[kv.first].push_back(nxi);
                nxi++;
            }
        }

        ex ut = 1;
        for(int i=0; i<Internal.nops(); i++) {
            ft = ft.expand();
            ft = subs(ft, Replacements, subs_options::algebraic);
            auto t2 = ft.coeff(Internal.op(i),2);
            auto t1 = ft.coeff(Internal.op(i),1);
            auto t0 = ft.subs(Internal.op(i)==0);
            ut *= t2;
            if(is_zero(t2)) return lst{0,0};
            ft = normal(t0-t1*t1/(4*t2));
        }
        ft = ex(0)-subs(ut*ft, Replacements, subs_options::algebraic);
        ut = subs(ut, Replacements, subs_options::algebraic);
        ex uf = normal(ut*ft);
        
        lst xRepl;
        for(int i=0; i<nxi; i++) xRepl.append(x(i));
        
        ex_is_less comp;
        for(auto pi : pgrp) {
            auto vi = pi.second;
            int nvi = vi.size();
            if(nvi<2) continue;
            ex uf1 = uf;
            auto xRepl1 = xRepl;
            Permutations(nvi, [comp,nvi,vi,uf,&uf1,xRepl,&xRepl1](const int *ns) {
                exmap x2x;
                for(int i=0; i<nvi; i++) x2x[x(vi[i])]=x(vi[ns[i]]);
                ex uf2 = uf.subs(x2x);
                if(comp(uf2, uf1)) {
                    uf1 = uf2;
                    xRepl1 = ex_to<lst>(subs(xRepl,x2x));
                }
            });
            uf = uf1;
            xRepl = xRepl1;
        }

        for(int i=0; i<nxi; i++) xRepl.let_op(i) = (x(i)==xRepl.op(i));
        ut = normal(ut.subs(xRepl));
        ft = normal(ft.subs(xRepl));
        return lst{ut, ft};
    }  
    
    /**
     * @brief Find Rules for Integrals or Master Integrals
     * @param fs vector of FIRE object
     * @param mi true for Master Integals
     * @return rules replacement and left integrals or left master integrals
     */
    pair<exmap,lst> FIRE::FindRules(vector<FIRE> & fs, bool mi) {
        exvector uf_mi_vec;
        if(mi) {
            uf_mi_vec = GiNaC_Parallel(-1, fs.size(), [mi,fs](int idx)->ex {
                const FIRE & fi = fs[idx]; // only here
                lst uf_mi_lst;
                for(auto mi : fi.MasterIntegrals) {
                    uf_mi_lst.append(lst{ fi.UF(mi.subs(F(w1,w2)==w2)), mi });
                }
                return uf_mi_lst;
            }, "MI");
        } else {
            uf_mi_vec = GiNaC_Parallel(-1, fs.size(), [mi,fs](int idx)->ex {
                const FIRE & fi = fs[idx]; // only here
                lst uf_mi_lst;
                for(auto mi : fi.Integrals) {
                    uf_mi_lst.append(lst{ fi.UF(mi), F(fi.ProblemNumber,mi) });
                }
                return uf_mi_lst;
            }, "I");
        }
    
        map<ex,lst,ex_is_less> group;
        int ntotal = 0;
        for(auto item1 : uf_mi_vec) {
            for(auto item : item1) {
                group[item.op(0)].append(item.op(1));
                ntotal++;
            }
        }

        exmap rules;
        lst int_lst;
        for(auto g : group) {
            lst gs = ex_to<lst>(g.second);
            for(int i=1; i<gs.nops(); i++) rules[gs.op(i)]=gs.op(0);
            int_lst.append(gs.op(0));
        }
        
        if(Verbose>2) cout << "  \\--FindRules: " << ntotal << " :> " << int_lst.nops() << " @ " << now(false) << endl;
        return make_pair(rules,int_lst);
    }
    
    /**
     * @brief Find Rules for Integrals or Master Integrals
     * @param fs vector of FIRE pointer object
     * @param mi true for Master Integals
     * @return rules replacement and left integrals or left master integrals
     */
    pair<exmap,lst> FIRE::FindRules(vector<FIRE*> & fs, bool mi) {
        exvector uf_mi_vec;
        if(mi) {
            uf_mi_vec = GiNaC_Parallel(-1, fs.size(), [mi,fs](int idx)->ex {
                const FIRE & fi = *(fs[idx]); // only here
                lst uf_mi_lst;
                for(auto mi : fi.MasterIntegrals) {
                    uf_mi_lst.append(lst{ fi.UF(mi.subs(F(w1,w2)==w2)), mi });
                }
                return uf_mi_lst;
            }, "MI");
        } else {
            uf_mi_vec = GiNaC_Parallel(-1, fs.size(), [mi,fs](int idx)->ex {
                const FIRE & fi = *(fs[idx]); // only here
                lst uf_mi_lst;
                for(auto mi : fi.Integrals) {
                    uf_mi_lst.append(lst{ fi.UF(mi), F(fi.ProblemNumber,mi) });
                }
                return uf_mi_lst;
            }, "II");
        }
    
        map<ex,lst,ex_is_less> group;
        int ntotal = 0;
        for(auto item1 : uf_mi_vec) {
            for(auto item : item1) {
                group[item.op(0)].append(item.op(1));
                ntotal++;
            }
        }

        exmap rules;
        lst int_lst;
        for(auto g : group) {
            lst gs = ex_to<lst>(g.second);
            for(int i=1; i<gs.nops(); i++) rules[gs.op(i)]=gs.op(0);
            int_lst.append(gs.op(0));
        }
        
        if(Verbose>2) cout << "  \\--FindRules: " << ntotal << " :> " << int_lst.nops() << " @ " << now(false) << endl;
        return make_pair(rules,int_lst);
    }

}
