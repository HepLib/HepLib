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

    static void a_print(const ex & ex_in, const print_context & c) {
        c.s << "a[" << ex_in << "]";
    }

    REGISTER_FUNCTION(a, do_not_evalf_params().print_func<print_dflt>(a_print))
    
    /**
     * @brief Do IBP Reduction 
     */
    void FIRE::Reduce() { 
        ProblemDimension = Propagators.nops();
        if(Integrals.nops()<1) return;
        int pn = 0; // to avoid unsigned short overflow in FIRE
        int pdim = ProblemDimension;
        lst InExternal;
        for(auto ii : Internal) InExternal.append(ii);
        for(auto ii : External) InExternal.append(ii);
        
        if(Pairs.nops()<1) {
            for(auto it : Internal) {
                for(auto ii : InExternal) Pairs.append(it*ii);
            }
            Pairs.sort();
            Pairs.unique();
            
            if(Pairs.nops()<pdim) {
                lst sps_ext;
                for(auto it : External) {
                    for(auto ii : External) sps_ext.append(it*ii);
                }
                sps_ext.sort();
                sps_ext.unique();
                for(auto item : sps_ext) {
                    auto item2 = subs_all(item,Replacements);
                    if(is_zero(item-item2)) Pairs.append(item);
                }
                Pairs.sort();
                Pairs.unique();
            }
        }
        
        if(Pairs.nops() != pdim) {
            cout << "Pairs = " << Pairs << endl;
            cout << "Propagators = " << Propagators << endl;
            throw Error("FIRE::Reduce: Pairs failed.");
        }
        
        lst sp2s, s2sp, ss;
        for(auto item : Pairs) {
            symbol si;
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst eqns, iWFs;
        for(int i=0; i<pdim; i++) {
            auto eq = Propagators.op(i).expand().subs(iEpsilon==0); // drop iEpsilon
            eq = eq.subs(sp2s, subs_options::algebraic);
            eq = eq.subs(Replacements, subs_options::algebraic);
            if(eq.has(iWF(w))) throw Error("iWF used in eq.");
            eqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(eqns, ss);
        if(s2p.nops() != pdim) throw Error("FIRE::Reduce: lsove failed.");

        IBPs.clear();
        exvector IBPvec;
        lst ns0;
        for(int i=0; i<pdim; i++) ns0.append(0);
        for(auto loop : Internal) {
            lst dp_lst;
            for(int i=0; i<pdim; i++) { 
                auto s = ex_to<Symbol>(loop);
                dp_lst.append(Propagators.op(i).diff(s));
            } 
            
            for(auto iep : InExternal) { 
                exmap nc_map;
                for(int i=0; i<pdim; i++) { // diff on each propagator
                    auto ns = ns0;
                    ns.let_op(i) = ns.op(i)+1; // note the covention
                    auto tmp = dp_lst.op(i) * iep;
                    tmp = mma_collect(tmp, InExternal);
                    tmp = tmp.subs(Replacements, subs_options::algebraic);
                    tmp = tmp.subs(sp2s, subs_options::algebraic);
                    tmp = tmp.subs(s2p, subs_options::algebraic);
                    tmp = ex(0) - a(i+1)*tmp;

                    for(int j=0; j<pdim; j++) {
                        auto cj = tmp.coeff(iWF(j));
                        if(is_zero(cj)) continue;
                        auto cns = ns;
                        cns.let_op(j) = cns.op(j)-1; // note the covention
                        nc_map[cns] = nc_map[cns] + cj;
                    }
                    tmp = tmp.subs(iWF(w)==0); // constant term
                    if(!is_zero(tmp)) nc_map[ns] = nc_map[ns] + tmp;
                }
                
                if(is_zero(loop-iep)) nc_map[ns0] = nc_map[ns0] + VectorDimension;
                bool ok = false;
                for(auto nc : nc_map) {
                    if(!is_zero(nc.second)) {
                        ok = true;
                        IBPvec.push_back(nc.second.subs(D==d)); // replace D to d
                    }
                }
                if(ok) IBPs.push_back(nc_map);
            }
        }

        auto Variables = gather_symbols(IBPvec);
        
        ostringstream start;

        start << "ExampleDimension[" << pn << "]=" << pdim << endl << endl;
        start << "ProblemNumber=" << pn << endl << endl;
        
        // .start - SBasisL
        if(Version==5) {
            PermutationsR(2, pdim, [pdim,pn,&start](const int *ns) {
                start << "SBasisL[" << pn << ",{";
                for(int i=0; i<pdim; i++) start << (ns[i]<1 ? -1 : 1) << (i<pdim-1 ? "," : "");
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
        for(int i=0; i<pdim; i++) start << (i+1) << (i<pdim-1 ? "," : "");
        start << "},{";
        // 1,1,1,...
        for(int i=0; i<pdim; i++) start << 1 << (i<pdim-1 ? "," : "");
        start << "},{";
        // 0,0,0,...
        for(int i=0; i<pdim; i++) start << 0 << (i<pdim-1 ? "," : "");
        start << "}}}" << endl << endl;
        
        // .start - SBasisR
        lst Rlst;
        Rlst.append(lst{});
        for(int i=0; i<pdim; i++) {
            let_op_append(Rlst, 0, -1);
        }
        for(auto lpi : Internal) {
            vector<int> ns_vec;
            lst ns0;
            for(int i=0; i<pdim; i++) ns0.append(1);
            for(int i=0; i<pdim; i++) {
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
            Rlst.remove_all();  // TODO: check this one
            for(auto cx : Cuts) {
                int ci = ex_to<numeric>(cx-1).to_int(); // start from 1 in Cuts
                lst ns0;
                for(int i=0; i<pdim; i++) ns0.append(1);
                ns0.let_op(ci) = -1;
                for(int n=0; n<std::pow(2,pdim-1); n++) {
                    int cn = n;
                    lst ns1 = ns0;
                    for(int j=0; j<pdim; j++) {
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
            for(int i=0; i<pdim; i++) start << iR.op(i) << (i<pdim-1 ? "," : "");
            start << "}]=True" << endl << endl;
        }
        
        // .start - Others
        start << "SBasisRL[" << pn << "]=0" << endl << endl;
        start << "HPI[" << pn << "]={}" << endl << endl;
        
        string sss = start.str();
        string_replace_all(sss, "=", " = ");
        string_replace_all(sss, ",", ", ");
        
        string spn = to_string(ProblemNumber);
        system(("mkdir -p " + WorkingDir).c_str());
        
        // .config
        ostringstream config;
        if(Version>5) config << "#compressor none" << endl;
        config << "#threads 1" << endl;
        config << "#fthreads 1" << endl;
        //config << "#fermat fer64" << endl;
        config << "#variables ";
        bool first = true;
        for(auto v : Variables) { 
            if(!islower(ex_to<symbol>(v).get_name()[0])) 
                throw Error("Reduce: Fermat requires a name must begin with a lower case letter.");
            config << (first ? "" : ",") << v; 
            first=false; 
        }
        config << endl;
        config << "#database db" << ProblemNumber << endl;
        if(Version>5 && pos_pref!=1) config << "#pos_pref "<< pos_pref << endl;
        if(Version>5) config << "#allIBP" << endl;
        if(Version==5) config << "#bucket 20" << endl;
        config << "#start" << endl;
        config << "#problem " << pn << " " << ProblemNumber << ".start" << endl;
        if(mi_pref.nops()>0) {
            ostringstream oss;
            oss << "{";
            int nn = mi_pref.nops();
            for(int i=0; i<nn; i++) oss << "{" << pn << "," << mi_pref.op(i) << (i<nn-1 ? "}," : "}");
            oss << "}";
            ofstream pref_out(WorkingDir+"/"+spn+".pref");
            pref_out << oss.str() << endl;
            pref_out.close();
            config << "#preferred " << ProblemNumber << ".pref" << endl;
        }
        config << "#integrals " << ProblemNumber << ".intg" << endl;
        config << "#output " << ProblemNumber << ".tables" << endl;
        
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
     
        // handle Cuts not equal 1, using #preferred
        if(Cuts.nops()>0 && mi_pref.nops()<1) {
            lst cIntegrals;
            auto mis = MasterIntegrals;
            MasterIntegrals.remove_all();
            for(auto item : mis) {
                lst mi = ex_to<lst>(item.op(1));
                bool isOK = true;
                for(auto cx : Cuts) {
                    if(!is_zero(mi.op(ex_to<numeric>(cx).to_int()-1)-1)) {
                        isOK = false;
                        break;
                    }
                }
                if(isOK) MasterIntegrals.append(item);
                else {
                    cIntegrals.append(mi);
                    auto pi = mi;
                    vector<int> ipos, ineg;
                    for(int i=0; i<mi.nops(); i++) {
                        bool isCut = false;
                        for(auto cx : Cuts) {
                            if(is_zero(cx-1-i)) {
                                isCut = true;
                                break;
                            }
                        }
                        if(isCut) pi.let_op(i)=1;
                        else if(mi.op(i)<=0) ineg.push_back(i);
                        else ipos.push_back(i);
                    }
                    
                    lst mi_pref_tmp;
                    int max = 2;
                    ex total = pow(numeric(max), ineg.size());
                    for(numeric in=0; in<total; in++) {
                        auto cin = in;
                        auto pi2 = pi;
                        for(int i=0; i<ineg.size(); i++) {
                            int re = mod(cin,max).to_int();
                            pi2.let_op(ineg[i]) = ex(0)-re;
                            cin = (cin-re)/max;
                        }
                        mi_pref_tmp.append(pi2);
                    }
                    
                    int max2 = 2*max+1;
                    total = pow(numeric(max2), ipos.size());
                    for(auto pi : mi_pref_tmp) {
                        for(numeric in=0; in<total; in++) {
                            auto cin = in;
                            auto pi2 = pi;
                            for(int i=0; i<ipos.size(); i++) {
                                int re = mod(cin,max2).to_int();
                                pi2.let_op(ipos[i]) = re-max;
                                cin = (cin-re)/max2;
                            }
                            mi_pref.append(pi2);
                        }
                    }
                }
            }
            
            // Reduce again
            mi_pref.sort();
            mi_pref.unique();
            if(mi_pref.nops()>0) {
                FIRE fire;
                fire.Propagators = Propagators;
                fire.Internal = Internal;
                fire.External = External;
                fire.Replacements = Replacements;
                fire.Pairs = Pairs;
                fire.ProblemNumber = ProblemNumber;
                fire.Cuts = Cuts;
                fire.Integrals = cIntegrals;
                fire.mi_pref = mi_pref;
                mi_pref.remove_all();
                fire.WorkingDir = WorkingDir + "_C"+to_string(ProblemNumber);
                fire.Reduce();
                system(("rm -rf " + fire.WorkingDir).c_str());
                for(auto item : fire.MasterIntegrals) MasterIntegrals.append(item);
                auto rules = Rules;
                Rules.remove_all();
                for(auto r : rules) Rules.append(r.op(0)==subs_naive(r.op(1),fire.Rules));
            }
        }
    }
    
    /**
     * @brief Sort for all permuations, and return xs w.r.t. 1st permutation
     * @param expr the input expression, as the sort key, no need of polynormial of xs
     * @param xs the permutation list
     * @return 1st of sorted xs
     */
    lst FIRE::SortPermutation(const ex & in_expr, const lst & xs) {
        auto expr = in_expr;
        bool isPoly = true;
        lst xRepl;
        map<ex,vector<int>,ex_is_less> pgrp;
        if(!expr.is_polynomial(xs)) {
            expr = expr.normal();
            expr = expr.numer_denom();
            if(!expr.is_polynomial(xs)) {
                isPoly = false;
                xRepl = xs;
                for(int i=0; i<xs.nops(); i++) pgrp[-1].push_back(i); // all permuations
            } else {
                expr = expr.op(0) * expr.op(1);
            }
        }
        
        if(isPoly) { // only for polynomials
            if(expr.has(coCF(w)) || expr.has(coVF(w))) throw Error("SortPermutation: coVF/coCF found.");
            expr = mma_collect(expr, xs, true, true);
            if(!is_a<add>(expr)) expr = lst{expr};
            exvector cvs;
            for(auto item : expr) {
                cvs.push_back(lst{
                    item.subs(lst{coCF(w)==w,coVF(w)==1}), // coefficient
                    item.subs(lst{coCF(w)==1,coVF(w)==w}) // xs monomial
                });
            }
            sort_vec_by(cvs,0);
                    
            int nxi = xs.nops();
            bool first = true;
            lst xkey[nxi];
            lst subkey[nxi];
            ex clast;
            for(auto cv : cvs) {
                ex cc = cv.op(0).expand();
                ex vv = cv.op(1);
                if(is_zero(cc)) continue;
                if(!first && !is_zero(cc-clast)) {
                    for(int i=0; i<nxi; i++) {
                        sort_lst(subkey[i]);
                        for(auto item : subkey[i]) xkey[i].append(item);
                        subkey[i].remove_all();
                    }
                } 
                first = false;
                clast = cc;
                for(int i=0; i<nxi; i++) subkey[i].append(vv.degree(xs.op(i)));
            }
            for(int i=0; i<nxi; i++) {
                sort_lst(subkey[i]);
                for(auto item : subkey[i]) xkey[i].append(item);
                subkey[i].remove_all();
            }
            
            exvector key_xi;
            for(int i=0; i<nxi; i++) {
                key_xi.push_back(lst{xkey[i], xs.op(i)});
                pgrp[xkey[i]].push_back(i); // i w.r.t. position of xs
            }
            sort_vec_by(key_xi,0);
            
            xRepl.remove_all();
            for(auto item : key_xi) xRepl.append(item.op(1));  
        }
        
        // pgrp - needs to permuation explicitly 
        expr = in_expr;
        for(auto pi : pgrp) {
            auto vi = pi.second;
            int nvi = vi.size();
            if(nvi<2) continue;
            ex expr1 = expr;
            auto xRepl1 = xRepl;
            Permutations(nvi, [xs,nvi,vi,expr,&expr1,xRepl,&xRepl1](const int *ns) {
                exmap x2x;
                for(int i=0; i<nvi; i++) x2x[xs.op(vi[i])]=xs.op(vi[ns[i]]);
                ex expr2 = subs_naive(expr,x2x);
                if(ex_less(expr2,expr1)) {
                    expr1 = expr2;
                    xRepl1 = ex_to<lst>(subs_naive(xRepl,x2x));
                }
            });
            expr = expr1;
            xRepl = xRepl1;
        }
        return xRepl;
    }
    
    /**
     * @brief UF function, from FIRE.m
     * @param idx exponent for the internal Propagator
     * @return lst of {U, F}
     */
    lst FIRE::LoopUF(const FIRE & fire, const ex & idx) {
        ex ft = 0;
        int nps = fire.Propagators.nops();
        int nxi=0;
        lst xs;
        for(int i=0; i<nps; i++) {
            if(is_zero(idx.op(i))) continue;
            if(!is_zero(idx.op(i)-1)) ft -= iWF(idx.op(i)) * x(nxi) * fire.Propagators.op(i);
            else ft -= x(nxi) * fire.Propagators.op(i);
            xs.append(x(nxi));
            nxi++;
        }
        
        ex ut = 1;
        for(int i=0; i<fire.Internal.nops(); i++) {
            ft = ft.expand();
            ft = subs_all(ft, fire.Replacements);
            auto t2 = ft.coeff(fire.Internal.op(i),2);
            auto t1 = ft.coeff(fire.Internal.op(i),1);
            auto t0 = ft.subs(fire.Internal.op(i)==0);
            ut *= t2;
            if(is_zero(t2)) return lst{0,0};
            ft = normal(t0-t1*t1/(4*t2));
        }
        ft = subs_all(ut*ft, fire.Replacements);
        ut = subs_all(ut, fire.Replacements);
        ex uf = normal(ut*ft);
        
        auto xRepl = SortPermutation(uf,xs);
        for(int i=0; i<nxi; i++) xRepl.let_op(i)=(xRepl.op(i)==x(i));
        ut = normal(subs_naive(ut,xRepl));
        ft = normal(subs_naive(ft,xRepl));
        return lst{ut, ft};
    }  
    
    /**
     * @brief UF function, from FIRE.m
     * @param idx exponent for the internal Propagator
     * @return lst of {U, F}
     */
    lst FIRE::UF(const ex & ps, const ex & ns, const ex & loops, const ex & tloops, const ex & lsubs, const ex & tsubs) {
        ex ft = 0;
        int nps = ps.nops();
        int nxi=0;
        lst xs;
        for(int i=0; i<nps; i++) {
            if(is_zero(ns.op(i))) continue;
            if(!is_zero(ns.op(i)-1)) ft -= iWF(ns.op(i)) * x(nxi) * ps.op(i);
            else ft -= x(nxi) * ps.op(i);
            xs.append(x(nxi));
            nxi++;
        }

        ex ut1 = 1;
        for(int i=0; i<loops.nops(); i++) {
            ft = ft.expand();
            ft = subs_all(ft, lsubs);
            auto t2 = ft.coeff(loops.op(i),2);
            auto t1 = ft.coeff(loops.op(i),1);
            auto t0 = ft.subs(loops.op(i)==0);
            ut1 *= t2;
            if(is_zero(t2)) return lst{0,0,0};
            ft = normal(t0-t1*t1/(4*t2));
        }
        ft = normal(subs_all(ut1*ft, lsubs));
        ut1 = normal(subs_all(ut1, lsubs));

        ex ut2 = 1;
        for(int i=0; i<tloops.nops(); i++) {
            ft = ft.expand();
            ft = subs_all(ft, tsubs);
            auto t2 = ft.coeff(tloops.op(i),2);
            auto t1 = ft.coeff(tloops.op(i),1);
            auto t0 = ft.subs(tloops.op(i)==0);
            ut2 *= t2;
            if(is_zero(t2)) return lst{0,0,0};
            ft = normal(t0-t1*t1/(4*t2));
        }
        ft = normal(subs_all(ut2*ft, tsubs));
        ut2 = normal(subs_all(ut2, tsubs));
        
        ex uf = normal(ut1*ut2*ft);
        
        lst xRepl = SortPermutation(uf,xs);
        for(int i=0; i<nxi; i++) xRepl.let_op(i)=(xRepl.op(i)==x(i));
        uf = uf.subs(xRepl);

        // z Permuatations
        if(tloops.nops()>1) {
            lst zs;
            auto nzi = tloops.nops();
            for(int i=0; i<nzi; i++) zs.append(z(i+1));
            auto zRepl = SortPermutation(uf,zs);
            for(int i=0; i<nzi; i++) zRepl.let_op(i)=(zRepl.op(i)==z(i+1));
            for(auto item : zRepl) xRepl.append(item);
        }

        ut1 = normal(ut1.subs(xRepl));
        ut2 = normal(ut2.subs(xRepl));
        ft = normal(ft.subs(xRepl));
        return lst{ut1, ut2, ft};
    }
    
    /**
     * @brief Find Rules for Integrals or Master Integrals
     * @param fs vector of FIRE object
     * @param mi true for Master Integals
     * @return rules replacement and left integrals or left master integrals
     */
    pair<exmap,lst> FIRE::FindRules(vector<FIRE> & fs, bool mi, std::function<lst(const FIRE &, const ex &)> uf) {
        exvector uf_mi_vec;
        if(mi) {
            uf_mi_vec = GiNaC_Parallel(fs.size(), [mi,fs,uf](int idx)->ex {
                const FIRE & fi = fs[idx]; // only here
                lst uf_mi_lst;
                for(auto mi : fi.MasterIntegrals) {
                    uf_mi_lst.append(lst{ uf(fi,mi.subs(F(w1,w2)==w2)), mi });
                }
                return uf_mi_lst;
            }, "MI");
        } else {
            uf_mi_vec = GiNaC_Parallel(fs.size(), [mi,fs,uf](int idx)->ex {
                const FIRE & fi = fs[idx]; // only here
                lst uf_mi_lst;
                for(auto mi : fi.Integrals) {
                    uf_mi_lst.append(lst{ uf(fi,mi), F(fi.ProblemNumber,mi) });
                }
                return uf_mi_lst;
            }, "FI");
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
    pair<exmap,lst> FIRE::FindRules(vector<FIRE*> & fs, bool mi, std::function<lst(const FIRE &, const ex &)> uf) {
        exvector uf_mi_vec;
        if(mi) {
            uf_mi_vec = GiNaC_Parallel(fs.size(), [mi,fs,uf](int idx)->ex {
                const FIRE & fi = *(fs[idx]); // only here
                lst uf_mi_lst;
                for(auto mi : fi.MasterIntegrals) {
                    uf_mi_lst.append(lst{ uf(fi,mi.subs(F(w1,w2)==w2)), mi });
                }
                return uf_mi_lst;
            }, "MI");
        } else {
            uf_mi_vec = GiNaC_Parallel(fs.size(), [mi,fs,uf](int idx)->ex {
                const FIRE & fi = *(fs[idx]); // only here
                lst uf_mi_lst;
                for(auto mi : fi.Integrals) {
                    uf_mi_lst.append(lst{ uf(fi,mi), F(fi.ProblemNumber,mi) });
                }
                return uf_mi_lst;
            }, "FI");
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
