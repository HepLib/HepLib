/**
 * @file
 * @brief IBP
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib {

    namespace {
        ex_is_less less;
    }

    static void a_print(const ex & ex_in, const print_context & c) {
        c.s << "a[" << ex_in << "]";
    }

    REGISTER_FUNCTION(a, do_not_evalf_params().print_func<print_dflt>(a_print))
    
    void IBP::Reduce() {
        Export();
        Run();
        Import();
    }
    
    void IBP::Export(string garfn) {
        archive ar;
        ar.archive_ex(Internal, "Internal");
        ar.archive_ex(External, "External");
        ar.archive_ex(Replacements, "Replacements");
        ar.archive_ex(Propagators, "Propagators");
        ar.archive_ex(Cuts, "Cuts");
        ar.archive_ex(DSP, "DSP");
        ar.archive_ex(ISP, "ISP");
        
        ar.archive_ex(Shift.size(), "NShift");
        int i = 0;
        for(auto kv : Shift) {
            ar.archive_ex(kv.first, ("ShiftK-"+to_string(i)).c_str());
            ar.archive_ex(kv.second, ("ShiftV-"+to_string(i)).c_str());
            i++;
        }
        
        ar.archive_ex(reCut, "reCut"); // bool
        ar.archive_ex(ProblemNumber, "ProblemNumber"); // int
        ar.archive_ex(Symbol(WorkingDir), "WorkingDir"); // string
        ar.archive_ex(PIntegrals, "PIntegrals");
        ar.archive_ex(MIntegrals, "MIntegrals");
        ar.archive_ex(Rules, "Rules");
        ar.archive_ex(IsAlwaysZero, "IsAlwaysZero"); // bool
        
        ofstream out(garfn);
        out << ar;
        out.close();
    }
    
    ex IBP::TO() {
        lst shift;
        for(auto kv : Shift) shift.append(lst{kv.first, kv.second});
        
        return lst{ Internal, External, Replacements, Propagators, Cuts, DSP, ISP, shift, reCut, ProblemNumber, Symbol(WorkingDir), PIntegrals, MIntegrals, Rules, IsAlwaysZero };
    }
    
    void IBP::Import(string garfn) {
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        
        Internal = ex_to<lst>(ar.unarchive_ex("Internal"));
        External = ex_to<lst>(ar.unarchive_ex("External"));
        Replacements = ex_to<lst>(ar.unarchive_ex("Replacements"));
        Propagators = ex_to<lst>(ar.unarchive_ex("Propagators"));
        Cuts = ex_to<lst>(ar.unarchive_ex("Cuts"));
        DSP = ex_to<lst>(ar.unarchive_ex("DSP"));
        ISP = ex_to<lst>(ar.unarchive_ex("ISP"));
        
        int n = ex_to<numeric>(ar.unarchive_ex("NShift")).to_int();
        for(int i=0; i<n; i++) {
            int key = ex_to<numeric>(ar.unarchive_ex(("ShiftK-"+to_string(i)).c_str())).to_int();
            ex val = ar.unarchive_ex(("ShiftV-"+to_string(i)).c_str());
            Shift[key] = val;
        }
        
        reCut = !(ar.unarchive_ex("reCut").is_zero()); // bool
        ProblemNumber = ex_to<numeric>(ar.unarchive_ex("ProblemNumber")).to_int(); // int
        WorkingDir = ex_to<Symbol>(ar.unarchive_ex("WorkingDir")).get_name();
        PIntegrals = ex_to<lst>(ar.unarchive_ex("PIntegrals"));
        MIntegrals = ex_to<lst>(ar.unarchive_ex("MIntegrals"));
        Rules = ex_to<lst>(ar.unarchive_ex("Rules"));
        IsAlwaysZero = !(ar.unarchive_ex("IsAlwaysZero").is_zero()); // bool
    
    }
    
    void IBP::FROM(ex s) {
        int i = 0;
        Internal = ex_to<lst>(s.op(i++));
        External = ex_to<lst>(s.op(i++));
        Replacements = ex_to<lst>(s.op(i++));
        Propagators = ex_to<lst>(s.op(i++));
        Cuts = ex_to<lst>(s.op(i++));
        DSP = ex_to<lst>(s.op(i++));
        ISP = ex_to<lst>(s.op(i++));
        lst shift = ex_to<lst>(s.op(i++));
        reCut = s.op(i++).is_equal(1);
        ProblemNumber = ex_to<numeric>(s.op(i++)).to_int();
        WorkingDir = ex_to<Symbol>(s.op(i++)).get_name();
        PIntegrals = ex_to<lst>(s.op(i++));
        MIntegrals = ex_to<lst>(s.op(i++));
        Rules = ex_to<lst>(s.op(i++));
        IsAlwaysZero = s.op(i++).is_equal(1);
        for(auto item : shift) Shift[ex_to<numeric>(item.op(0)).to_int()] = item.op(1);
    }
    
    void IBP::ReShare(const vector<IBP*> & fs) {
        string garfn = to_string(getpid())+".ReShare.gar";
        if(true) {
            archive ar;
            for(int i=0; i<fs.size(); i++) {
                string si = to_string(i);
                ar.archive_ex(fs[i]->Internal, (si+"-1").c_str());
                ar.archive_ex(fs[i]->External, (si+"-2").c_str());
                ar.archive_ex(fs[i]->Replacements, (si+"-3").c_str());
                ar.archive_ex(fs[i]->Propagators, (si+"-4").c_str());
                ar.archive_ex(fs[i]->DSP, (si+"-5").c_str());
                ar.archive_ex(fs[i]->ISP, (si+"-6").c_str());
                ar.archive_ex(fs[i]->PIntegrals, (si+"-7").c_str());
                ar.archive_ex(fs[i]->MIntegrals, (si+"-8").c_str());
                ar.archive_ex(fs[i]->Rules, (si+"-9").c_str());
            }
            ofstream out(garfn);
            out << ar;
            out.close();
        }
        if(true) {
            archive ar;
            ifstream in(garfn);
            in >> ar;
            in.close();
            
            map<string, ex> dict;
            for(int i=0; i<ar.num_expressions(); i++) {
                string name;
                ex res = ar.unarchive_ex(name, i);
                dict[name] = res;
            }
            
            for(int i=0; i<fs.size(); i++) {
                string si = to_string(i);
                fs[i]->Internal = ex_to<lst>(dict[si+"-1"]);
                fs[i]->External = ex_to<lst>(dict[si+"-2"]);
                fs[i]->Replacements = ex_to<lst>(dict[si+"-3"]);
                fs[i]->Propagators = ex_to<lst>(dict[si+"-4"]);
                fs[i]->DSP = ex_to<lst>(dict[si+"-5"]);
                fs[i]->ISP = ex_to<lst>(dict[si+"-6"]);
                fs[i]->PIntegrals = ex_to<lst>(dict[si+"-7"]);
                fs[i]->MIntegrals = ex_to<lst>(dict[si+"-8"]);
                fs[i]->Rules = ex_to<lst>(dict[si+"-9"]);
            }
        }
    }
    
    pair<exmap,lst> IBP::FindRules(bool is_mi) {
        vector<IBP*> ibps;
        ibps.push_back(this);
        auto rs_mis = HepLib::FindRules(ibps, is_mi);
        if(is_mi && rs_mis.first.size()>0) {
            auto nr = Rules.nops();
            for(int i=0; i<nr; i++) {
                auto ri = Rules.op(i);
                Rules.let_op(i) = (ri.op(0) == ri.op(1).subs(rs_mis.first, nopat));
            }
            for(auto mi : MIntegrals) {
                auto mi2 = mi.subs(rs_mis.first, nopat);
                if(mi==mi2) continue;
                Rules.append(mi==mi2);
            }
            MIntegrals = rs_mis.second;
            sort_lst(MIntegrals);
        }
        return rs_mis;
    }
        
    /**
     * @brief Sort for all permuations, and return xs w.r.t. 1st permutation
     * @param in_expr the input expression, as the sort key, no need of polynormial of xs
     * @param xs the permutation list
     * @return x-replacement
     */
    exmap SortPermutation(const ex & in_expr, const lst & xs) {
        auto expr = in_expr;
        bool isPoly = true;
        exmap xmap; // note that we need xmap.size == xs.nops()
        
        map<ex,vector<int>,ex_is_less> pgrp;
        if(!expr.is_polynomial(xs)) {
            if(Verbose>2) cout << "SortPermutation: NOT a polynomials, try ALL permutations!" << endl;
            expr = expr.numer_denom();
            if(true || !expr.op(0).is_polynomial(xs) || !expr.op(1).is_polynomial(xs)) {
                isPoly = false;
                for(int i=0; i<xs.nops(); i++) {
                    pgrp[-1].push_back(i); // all permuations
                    xmap[xs.op(i)] = xs.op(i); 
                }
            } else {
                expr = expr.op(0) * expr.op(1);
            }
        }
        
        if(isPoly) { // only for polynomials
            auto cv_lst = collect_lst(expr, xs);
            exvector cvs;
            for(auto item : cv_lst) cvs.push_back(item);
            sort_vec(cvs); // first sort by coefficient
                    
            int nxi = xs.nops();
            bool first = true;
            lst xkey[nxi];
            lst subkey[nxi];
            ex clast;
            for(auto cv : cvs) {
                ex cc = cv.op(0);
                ex vv = cv.op(1);
                if(is_zero(cc)) continue;
                if(!first && !is_zero(cc-clast)) { // equal coefficient to the same key
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

            for(int i=0; i<nxi; i++) { // add last subkey
                sort_lst(subkey[i]); // need sort due to equal coefficient
                for(auto item : subkey[i]) xkey[i].append(item);
                subkey[i].remove_all();
            }
            
            exvector key_xi;
            for(int i=0; i<nxi; i++) key_xi.push_back(lst{xkey[i], xs.op(i)});
            sort_vec(key_xi); // first sort by key
            
            for(int i=0; i<nxi; i++) {
                auto xi = key_xi[i].op(1);
                auto xki = key_xi[i].op(0);
                xmap[xi] = xs.op(i); // initial x-replacement
                pgrp[xki].push_back(i);
            }
            expr = expr.subs(xmap, nopat); // need to update expr before permutations
            expr = collect_ex(expr, xs); // need collect here
        }
 
        // pgrp - needs to permuation explicitly 
        ex expr_min = expr;
        exmap xmap_min;
        long long npt = 1;
        for(auto pi : pgrp) {
            int nvi = pi.second.size();
            for(int i=1; i<=nvi; i++) npt *= i;  // nvi!
        }
        for(long long np=0; np<npt; np++) {
            long long npc = np;
            exmap xmap_tmp;
            for(auto pi : pgrp) {
                auto vi = pi.second;
                int nvi = vi.size();
                if(nvi<2) continue;
                long long npti = 1;
                for(int i=1; i<=nvi; i++) npti *= i; // nvi!
                int ck = npc % npti; // kth permutation
                npc = npc / npti;
            
                // https://stackoverflow.com/questions/1995328/are-there-any-better-methods-to-do-permutation-of-string
                auto vip = vi;
                int k = ck;
                for(int j=1; j<nvi; ++j) { // startoverflow: j starts from 1
                    std::swap(vip[k%(j+1)], vip[j]); 
                    k=k/(j+1);
                }
                
                for(int j=0; j<nvi; j++) xmap_tmp[xs.op(vi[j])] = xs.op(vip[j]);
            }
            ex expr_tmp = expr.subs(xmap_tmp,nopat);
            if(ex_less(expr_tmp, expr_min)) {
                expr_min = expr_tmp;
                xmap_min = xmap_tmp;
            }
        }

        for(auto & kv : xmap) kv.second = kv.second.subs(xmap_min,nopat);
        return xmap;
    }
    
    /**
     * @brief UF function
     * @param IBP the IBP object
     * @param idx exponent for the internal Propagator
     * @return lst of {U, F, sign}
     */
    lst LoopUF(const IBP & ibp, const ex & idx) {
        
        auto props = ibp.Propagators;
        
        // handle sign
        ex sign = 1;
        for(int i=0; i<props.nops(); i++) {
            auto ipr = expand_ex(props.op(i), ibp.Internal);
            
            if(ipr.has(iEpsilon)) {
                auto cc = ipr.coeff(iEpsilon);
                if(is_a<numeric>(cc)) {
                    if(cc<0) { // using +iEpsilon
                        sign *= pow(-1, idx.op(i));
                        props.let_op(i) = ex(0)-props.op(i);
                    }
                    props.let_op(i) = props.op(i).subs(iEpsilon==0,nopat);
                    goto sign_done;
                } else throw Error("LoopUF: sign of iEpsilon NOT determined.");
            }
            
            for(auto lp : ibp.Internal) {
                if(ipr.degree(lp)==2) {
                    auto cc = ipr.coeff(lp,2);
                    if(is_a<numeric>(cc)) {
                        if(cc<0) { // using +l^2
                            sign *= pow(-1, idx.op(i));
                            props.let_op(i) = ex(0)-props.op(i);
                        }
                        goto sign_done;
                    }
                }
            }
            
            sign_done: ;
        }
        
        ex ut, ft, uf;
        lst key;
        lst xs;
        exmap x2ax;
        int nxi=0;
        int nps = props.nops();
        ft = 0;
        for(int i=0; i<nps; i++) {
            if(is_zero(idx.op(i))) {
                key.append(0);
                continue;
            }
            key.append(1);
            if(!is_zero(idx.op(i)-1)) x2ax[x(nxi)] = a(idx.op(i)) * x(nxi);
            xs.append(x(nxi));
            ft -= x(nxi) * props.op(i); // only used when no cache item
            nxi++;
        }
        
        static map<ex,exmap,ex_is_less> cache_by_prop;
        if(using_cache && cache_limit>0 && cache_by_prop.size() > cache_limit/10) cache_by_prop.clear();
        exmap & cache = cache_by_prop[lst{props,ibp.Internal}];
        if(!using_cache || cache.find(key)==cache.end()) { // no cache item
            ut = 1;
            ft = expand_ex(ft);
            ft = subs_all(ft, ibp.Replacements);
            for(int i=0; i<ibp.Internal.nops(); i++) {
                auto t2 = ft.coeff(ibp.Internal.op(i),2);
                auto t1 = ft.coeff(ibp.Internal.op(i),1);
                auto t0 = ft.subs(ibp.Internal.op(i)==0,nopat);
                ut *= t2;
                if(is_zero(t2)) return lst{0,0,1};
                ft = exnormal(t0-t1*t1/(4*t2));
                ft = expand_ex(ft);
                ft = subs_all(ft, ibp.Replacements);
            }
            ft = exnormal(ut*ft);
            ft = exnormal(subs_all(ft, ibp.Replacements));
            ut = exnormal(subs_all(ut, ibp.Replacements));
            uf = exnormal(ut+ft); // ut*ft, replay with ut+ft
            
            if(using_cache) cache[key] = lst{ut,ft,uf};
        } else {
            auto cc = cache[key];
            ut = cc.op(0);
            ft = cc.op(1);
            uf = cc.op(2);
        }

        ut = ut.subs(x2ax,nopat);
        ft = ft.subs(x2ax,nopat);
        uf = uf.subs(x2ax,nopat);
        
        uf = uf.subs(MapPreSP);
        auto xmap = SortPermutation(uf,xs);
        ut = (ut.subs(xmap,nopat));
        ft = (ft.subs(xmap,nopat));
        return lst{ut, ft, sign};
    }  
    
    /**
     * @brief UF function, from FIRE.m
     * @param ps the list of propagator
     * @param ns the list of exponent
     * @param loops the list of loop momenta
     * @param tloops the list of transverse/quasi momenta
     * @param lsubs the replacements for loops
     * @param tsubs the replacements for tloops
     * @return lst of {U1, U2, F, sign}
     */
    lst UF(const ex & props, const ex & ns, const ex & loops, const ex & tloops, const ex & lsubs, const ex & tsubs) {
        auto ps = props;
        
        // handle sign
        ex sign = 1;
        for(int i=0; i<ps.nops(); i++) {
            auto ipr = expand_ex(ps.op(i), loops);
            
            if(ipr.has(iEpsilon)) {
                auto cc = ipr.coeff(iEpsilon);
                if(is_a<numeric>(cc)) {
                    if(cc<0) { // using +iEpsilon
                        sign *= pow(-1, ns.op(i));
                        ps.let_op(i) = ex(0)-ps.op(i);
                    }
                    ps.let_op(i) = ps.op(i).subs(iEpsilon==0,nopat);
                    goto sign_done;
                } else throw Error("UF: sign of iEpsilon NOT determined.");
            }
            
            for(auto lp : loops) {
                if(ipr.degree(lp)==2) {
                    auto cc = ipr.coeff(lp,2);
                    if(is_a<numeric>(cc)) {
                        if(cc<0) { // using +l^2
                            sign *= pow(-1, ns.op(i));
                            ps.let_op(i) = ex(0)-ps.op(i);
                        }
                        goto sign_done;
                    }
                }
            }
            
            ipr = expand(ipr).subs(lsubs).subs(tsubs);
            for(auto lp : tloops) {
                if(ipr.degree(lp)==2) {
                    auto cc = ipr.coeff(lp,2);
                    if(is_a<numeric>(cc)) {
                        if(cc<0) { // using +l^2
                            sign *= pow(-1, ns.op(i));
                            ps.let_op(i) = ex(0)-ps.op(i);
                        }
                        goto sign_done;
                    }
                }
            }
            
            sign_done: ;
        }
        
        ex ut1, ut2, ft, uf;
        lst key;
        lst xs;
        exmap x2ax;
        int nxi=0;
        int nps = ps.nops();
        ft = 0;
        for(int i=0; i<nps; i++) {
            if(is_zero(ns.op(i))) {
                key.append(0);
                continue;
            }
            key.append(1);
            if(!is_zero(ns.op(i)-1)) x2ax[x(nxi)] = a(ns.op(i)) * x(nxi);
            xs.append(x(nxi));
            ft -= x(nxi) * ps.op(i); // only used when no cache item
            nxi++;
        }
        
        static map<ex,exmap,ex_is_less> cache_by_prop;
        exmap & cache = cache_by_prop[lst{ps,loops,tloops}];
        if(!using_cache || cache.find(key)==cache.end()) { // no cache item
            ut1 = 1;
            ft = expand(ft);
            ft = subs_all(ft, lsubs);
            for(int i=0; i<loops.nops(); i++) {
                auto t2 = ft.coeff(loops.op(i),2);
                auto t1 = ft.coeff(loops.op(i),1);
                auto t0 = ft.subs(loops.op(i)==0,nopat);
                ut1 *= t2;
                if(is_zero(t2)) return lst{0,0,0,1};
                ft = exnormal(t0-t1*t1/(4*t2));
                ft = expand(ft);
                ft = subs_all(ft, lsubs);
            }
            ft = exnormal(ut1*ft);
            ft = exnormal(subs_all(ft, lsubs));
            ut1 = exnormal(subs_all(ut1, lsubs));

            ut2 = 1;
            ft = expand(ft);
            ft = subs_all(ft, tsubs);
            for(int i=0; i<tloops.nops(); i++) {
                auto t2 = ft.coeff(tloops.op(i),2);
                auto t1 = ft.coeff(tloops.op(i),1);
                auto t0 = ft.subs(tloops.op(i)==0,nopat);
                ut2 *= t2;
                if(is_zero(t2)) return lst{0,0,0,1};
                ft = exnormal(t0-t1*t1/(4*t2));
                ft = expand(ft);
                ft = subs_all(ft, tsubs);
            }
            ft = exnormal(ut2*ft);
            ft = exnormal(subs_all(ft, tsubs));
            ut2 = exnormal(subs_all(ut2, tsubs));
            
            uf = exnormal(ut1*ut2*ft);
            if(using_cache) cache[key] = lst{ut1,ut2,ft,uf};
        } else {
            auto cc = cache[key];
            ut1 = cc.op(0);
            ut2 = cc.op(1);
            ft = cc.op(2);
            uf = cc.op(3);
        }
        ut1 = ut1.subs(x2ax,nopat);
        ut2 = ut2.subs(x2ax,nopat);
        ft = ft.subs(x2ax,nopat);
        uf = uf.subs(x2ax,nopat);
        
        uf = uf.subs(MapPreSP);
        auto xmap = SortPermutation(uf,xs);
        uf = uf.subs(xmap,nopat);

        // z Permuatations
        if(tloops.nops()>1) {
            lst zs;
            auto nzi = tloops.nops();
            for(int i=0; i<nzi; i++) zs.append(z(i+1));
            auto zmap = SortPermutation(uf,zs);
            for(auto kv : zmap) xmap[kv.first] = kv.second; // add to xmap
        }

        ut1 = (ut1.subs(xmap,nopat));
        ut2 = (ut2.subs(xmap,nopat));
        ft = (ft.subs(xmap,nopat));
        return lst{ut1, ut2, ft, sign};
    }
    
    /**
     * @brief Find Rules for Integrals or Master Integrals
     * @param fs vector of IBP pointer object
     * @param mi true for Master Integals
     * @param uf the function to compute the UF polynomial
     * @return rules replacement and left integrals or left master integrals
     */
    pair<exmap,lst> FindRules(vector<IBP*> fs, bool mi, std::function<lst(const IBP &, const ex &)> uf) {
        vector<pair<IBP*,ex>> ibp_idx_vec;
        if(mi) for(auto fi : fs) for(auto item : fi->MIntegrals) ibp_idx_vec.push_back(make_pair(fi, item));
        else for(auto fi : fs) for(auto item : fi->Integrals) ibp_idx_vec.push_back(make_pair(fi, F(fi->ProblemNumber,item)));
        
        exvector uf_smi_vec = GiNaC_Parallel(ibp_idx_vec.size(), [&ibp_idx_vec,&uf](int idx)->ex {
            auto p = ibp_idx_vec[idx];
            const IBP & fi = (*p.first);
            auto mi = p.second;
            auto ks = uf(fi,mi.subs(F(w1,w2)==w2));
            int nk = ks.nops()-1;
            lst key;
            for(int i=0; i<nk; i++) key.append(expand(ks.op(i)));
            lst pi = fi.Propagators; 
            return lst{ key, lst{ pi, ks.op(nk), mi } };  // ks.op(nk) -> the sign
        }, "FR");
            
        map<ex,lst,ex_is_less> group;
        int ntotal = 0;
        for(auto item : uf_smi_vec) {
            group[item.op(0)].append(item.op(1));
            ntotal++;
        }

        exmap rules;
        lst int_lst;
        exset pis;
        if(true) { // single element case
            exset ks, vs;
            for(auto g : group) {
                if(g.second.nops()==1) {
                    ks.insert(g.first);
                    vs.insert(g.second.op(0));
                }
            }
            for(auto vi : vs) {
                pis.insert(vi.op(0));
                int_lst.append(vi.op(2));
            }
            for(auto ki : ks) group.erase(ki);
        }
        
        while(!group.empty()) {
            if(pis.size()>0) {
                exset ks;
                for(auto g : group) {
                    ex c0, v0;
                    for(auto gi : g.second) {
                        for(auto pi : pis) {
                            if(gi.op(0).is_equal(pi)) {
                                ks.insert(g.first);
                                c0 = gi.op(1);
                                v0 = gi.op(2);
                                goto found;
                            }
                        }
                    }
                    continue; // if not found pi
                    found: ;
                    int_lst.append(v0);
                    for(auto gi : g.second) {
                        auto ci = gi.op(1);
                        auto vi = gi.op(2);
                        if(v0.is_equal(vi)) continue;
                        rules[vi] = v0 * c0 / ci;
                    }
                }
                for(auto ki : ks) group.erase(ki);
                if(group.empty()) break;
            }
            pis.clear();
            
            int cur = 0;
            for(auto g : group) {
                cur++;
                if(cur>100) break;
                pis.insert(g.second.op(0).op(0));
            }
        }
        
        if(Verbose>2) cout << "  \\--FindRules: " << WHITE << ntotal << " :> " << int_lst.nops() << RESET << " @ " << now(false) << endl;
        return make_pair(rules,int_lst);
    }
    
    static matrix Redrowech(const matrix & mat) {
        static map<pid_t, Fermat> fermat_map;
        static int v_max = 0;

        auto pid = getpid();
        if((fermat_map.find(pid)==fermat_map.end())) { // init section
            fermat_map[pid].Init();
            v_max = 0;
        }
        Fermat &fermat = fermat_map[pid];
        
        lst rep_vs;
        ex tree = mat;
        for(const_preorder_iterator i = tree.preorder_begin(); i != tree.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e) || e.match(a(w))) rep_vs.append(e);
        }
        rep_vs.sort();
        rep_vs.unique();
        sort_lst(rep_vs);
        
        exmap v2f;
        symtab st;
        int fvi = 0;
        for(auto vi : rep_vs) {
            auto name = "v" + to_string(fvi);
            v2f[vi] = Symbol(name);
            st[name] = vi;
            fvi++;
        }
        
        stringstream ss;
        if(fvi>111) {
            cout << rep_vs << endl;
            throw Error("Fermat: Too many variables.");
        }
        if(fvi>v_max) {
            for(int i=v_max; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            v_max = fvi;
        }
        
        int nrow = mat.rows();
        int ncol = mat.cols();
        
        ss << "Array m[" << nrow << "," << ncol << "];" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        
        ss << "[m]:=[(";
        for(int c=0; c<ncol; c++) {
            for(int r=0; r<nrow; r++) {
                ss << mat(r,c).subs(iEpsilon==0,nopat).subs(v2f,nopat) << ",";
            }
        }
        ss << ")];" << endl;
        ss << "Redrowech([m]);" << endl;
        auto tmp = ss.str();
        string_replace_all(tmp,",)]",")]");
        fermat.Execute(tmp);
        ss.clear();
        ss.str("");

        ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
        ss << "![m" << endl;
        auto ostr = fermat.Execute(ss.str());
        ss.clear();
        ss.str("");
        //fermat.Exit();
        
        // note the order, before exfactor (normal_fermat will be called again here)
        ss << "&(U=0);" << endl; // disable ugly printing
        ss << "@([m]);" << endl;
        ss << "&_G;" << endl;
        fermat.Execute(ss.str());
        ss.clear();
        ss.str("");

        // make sure last char is 0
        if(ostr[ostr.length()-1]!='0') throw Error("Direc::Export, last char is NOT 0.");
        ostr = ostr.substr(0, ostr.length()-1);
        string_trim(ostr);
        
        ostr.erase(0, ostr.find(":=")+2);
        string_replace_all(ostr, "[", "{");
        string_replace_all(ostr, "]", "}");
        Parser fp(st);
        matrix mr(nrow, ncol);
        auto res = fp.Read(ostr);
        for(int r=0; r<nrow; r++) {
            auto cur = res.op(r);
            for(int c=0; c<ncol; c++) mr(r,c) = cur.op(c);
        }
        return mr;
    }
    
    bool IBP::IsZero(ex sector) {
        auto props = Propagators;
        lst xs;
        int nxi=0;
        int nps = props.nops();
        ex ft = 0;
        for(int i=0; i<nps; i++) {
            if(is_zero(sector.op(i))) continue;
            xs.append(x(nxi));
            ft -= x(nxi) * props.op(i);
            nxi++;
        }
        
        ex ut = 1;
        ft = expand(ft);
        ft = subs_all(ft, Replacements);
        for(int i=0; i<Internal.nops(); i++) {
            auto t2 = ft.coeff(Internal.op(i),2);
            auto t1 = ft.coeff(Internal.op(i),1);
            auto t0 = ft.subs(Internal.op(i)==0,nopat);
            ut *= t2;
            if(is_zero(t2)) return true;
            ft = normal(t0-t1*t1/(4*t2));
            ft = expand(ft);
            ft = subs_all(ft, Replacements);
        }
        ut = normal(subs_all(ut, Replacements));
        ft = normal(ut*ft);
        ft = normal(subs_all(ft, Replacements));
    
        ex G = ut + ft;
        ex sum = 0;
        lst ks;
        exmap ks20;
        for(auto xi : xs) {
            symbol ki;
            ks.append(ki);
            ks20[ki] = 0;
            sum += ki * xi * diff_ex(G,xi);
        }
        sum -= G;
        auto cvs = collect_lst(sum,x(w));
        int rn = cvs.nops();
        int cn = ks.nops();
        matrix mat(rn,cn+1);
        int ri = 0;
        for(auto cv : cvs) {
            int ci = 0;
            for(auto ki : ks) {
                mat(ri,ci) = cv.op(0).coeff(ki);
                ci++;
            }
            mat(ri,cn) = cv.op(0).subs(ks20,nopat);
            ri++;
        }
        auto mat2 = Redrowech(mat);
        for(int ri=rn-1; ri>=0; ri--) {
            for(int ci=0; ci<cn; ci++) {
                if(!is_zero(mat2(ri,ci))) return true;
            }
            if(!is_zero(mat2(ri,cn))) return false;
        }
        return true;
    }
    
    exmap IBP::SP2Pn() { // sp -> {c0,c1,c2,...,cn} coefficient of each propagator, c0 is remaining constant
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
        
        int pdim = Propagators.nops();
        if(ISP.nops() > pdim) {
            cout << "ISP = " << ISP << endl;
            cout << "Propagators = " << Propagators << endl;
            throw Error("IBP::SP2Pn: #(ISP) > #(Propagators).");
        }
        
        lst sp2s, s2sp, ss;
        int _pic=0;
        for(auto item : ISP) {
            _pic++;
            symbol si("sp"+to_string(_pic));
            ss.append(si);
            sp2s.append(item==si);
            s2sp.append(si==item);
        }
        
        lst eqns;
        for(int i=0; i<ISP.nops(); i++) { // note NOT pdim
            auto eq = expand(Propagators.op(i)).subs(iEpsilon==0); // drop iEpsilon
            eq = eq.subs(sp2s, algbr);
            eq = eq.subs(Replacements, algbr);
            if(eq.has(iWF(w))) throw Error("IBP::SP2Pn, iWF used in eq.");
            eqns.append(eq == iWF(i));
        }
        auto s2p = lsolve(eqns, ss);
        if(s2p.nops() != ISP.nops()) {
            cout << ISP << endl << s2p << endl << eqns << endl;
            throw Error("IBP::SP2Pn: lsolve failed.");
        }
    
        exmap smap;
        for(auto r : s2p) {
            ex k = r.op(0).subs(s2sp, nopat);
            ex v = r.op(1);
            ex chk = v.subs(iWF(w)==0);
            lst res;
            res.append(chk);
            for(int i=0; i<ISP.nops(); i++) {
                ex t = v.coeff(iWF(i));
                res.append(t);
                chk += t*iWF(i);
            }
            if(!normal(chk-v).is_zero()) throw Error("IBP::SP2Pn check failed.");
            smap[k] = res;
        }
        return smap;
    }
    
    exmap IBP::Dinv(const lst & ns) {
        lst InExternal;
        for(auto ii : Internal) InExternal.append(ii);
        for(auto ii : External) InExternal.append(ii);
        auto & es = External;
        int eN = es.nops();
        int pN = Propagators.nops();
        matrix G(eN, eN);
        for(int r=0; r<eN; r++) for(int c=0; c<eN; c++) G(r,c) = (es.op(r)*es.op(c)).subs(Replacements,algbr);
        matrix Gi = G.inverse();
        // partial J/partial pi^2
        exmap spmap;
        auto sp2pn = SP2Pn();
        for(int p1i=0; p1i<eN; p1i++) {
        for(int p2i=p1i; p2i<eN; p2i++) {
            ex p1 = es.op(p1i);
            ex p2 = es.op(p2i);
            ex pf = 1;
            if(p1i==p2i) pf = 1/ex(2);
            ex res = 0;
            for(int pi=0; pi<pN; pi++) {
                lst ns2 = ns;
                ns2.let_op(pi) = ns.op(pi)+1;
                ex dpi = -ns.op(pi)*diff_ex(Propagators.op(pi), p1);
                for(int i=0; i<eN; i++) {
                    ex idpi = expand_ex(dpi*es.op(i)).subs(Replacements,algbr);
                    auto cvs = collect_lst(idpi, InExternal);
                    for(auto cv : cvs) {
                        if(is_zero(cv.op(1)-1)) res += cv.op(0)*pf*Gi(i,p2i)*F(ProblemNumber, ns2);
                        else {
                            auto f = sp2pn.find(cv.op(1));
                            if(f==sp2pn.end()) throw Error("IBP::DExt, Not found.");
                            auto cps = f->second;
                            res += cv.op(0)*pf*Gi(i,p2i)*cps.op(0)*F(ProblemNumber, ns2);
                            for(int j=1; j<pN+1; j++) {
                                if(is_zero(cps.op(j))) continue;
                                lst ns3 = ns2;
                                ns3.let_op(j-1) = ns2.op(j-1)-1;
                                res += cv.op(0)*pf*Gi(i,p2i)*cps.op(j)*F(ProblemNumber, ns3);
                            }
                        }
                    }
                }
            }
            spmap[p1*p2] = collect_ex(res,F(w1,w2));
        }}
        return spmap;
    }
    
    ex IBP::D(const ex & x, const lst & ns) {
        lst InExternal;
        for(auto ii : Internal) InExternal.append(ii);
        for(auto ii : External) InExternal.append(ii);
        int pN = Propagators.nops();
        auto sp2pn = SP2Pn();
        
        ex res = 0;
        for(int pi=0; pi<pN; pi++) { // Direct Diff for each Propagator
            lst ns2 = ns;
            ns2.let_op(pi) = ns.op(pi)+1;
            ex Pi = Propagators.op(pi);
            Pi = Pi.subs(Replacements);
            ex dpi = -ns.op(pi)*diff_ex(Pi, x);
            auto cvs = collect_lst(dpi, InExternal);
            for(auto cv : cvs) {
                if(is_zero(cv.op(1)-1)) res += cv.op(0)*F(ProblemNumber, ns2);
                else {
                    auto f = sp2pn.find(cv.op(1));
                    if(f==sp2pn.end()) throw Error("IBP::DExt, Not found.");
                    auto cps = f->second;
                    res += cv.op(0)*cps.op(0)*F(ProblemNumber, ns2);
                    for(int j=1; j<pN+1; j++) {
                        if(is_zero(cps.op(j))) continue;
                        lst ns3 = ns2;
                        ns3.let_op(j-1) = ns2.op(j-1)-1;
                        res += cv.op(0)*cps.op(j)*F(ProblemNumber, ns3);
                    }
                }
            }
        }

        // InDirect Diff from external SP
        auto dsp = Dinv(ns);
        auto eN = External.nops();
        for(int i=0; i<eN; i++) {
        for(int j=i; j<eN; j++) {
            ex sp = External.op(i) * External.op(j);
            auto f = dsp.find(sp);
            if(f==dsp.end()) throw Error("DESS::InitDE, sp NOT found.");
            auto rsp = sp.subs(Replacements, algbr);
            if(sp==rsp) throw Error("DESS::InitDE, sp==rsp, Replacements NOT work.");
            res += f->second * diff_ex(rsp, x);
        }}

        return res;
    }
    
    void IBP::RM(bool keep_start_config) {
        string spn = to_string(ProblemNumber);
        if(!keep_start_config) {
            file_remove(WorkingDir+"/"+spn+".start");
            file_remove(WorkingDir+"/"+spn+".config");
        }
        file_remove(WorkingDir+"/"+spn+".intg");
        file_remove(WorkingDir+"/"+spn+".tables");
    }
    
    ex GPolynomial(const IBP & ibp) {
        
        auto props = ibp.Propagators;
        ex ut, ft, uf;
        lst key;
        lst xs;
        exmap x2ax;
        int nxi=0;
        int nps = props.nops();
        ft = 0;
        for(int i=0; i<nps; i++) {
            xs.append(x(nxi));
            ft -= x(nxi) * props.op(i); // only used when no cache item
            nxi++;
        }
        
        ut = 1;
        ft = expand_ex(ft);
        ft = subs_all(ft, ibp.Replacements);
        for(int i=0; i<ibp.Internal.nops(); i++) {
            auto t2 = ft.coeff(ibp.Internal.op(i),2);
            auto t1 = ft.coeff(ibp.Internal.op(i),1);
            auto t0 = ft.subs(ibp.Internal.op(i)==0,nopat);
            ut *= t2;
            if(is_zero(t2)) return lst{0,0,1};
            ft = exnormal(t0-t1*t1/(4*t2));
            ft = expand_ex(ft);
            ft = subs_all(ft, ibp.Replacements);
        }
        ft = exnormal(ut*ft);
        ft = exnormal(subs_all(ft, ibp.Replacements));
        ut = exnormal(subs_all(ut, ibp.Replacements));
        uf = exnormal(ut+ft); // ut*ft, replay with ut+ft
            
        uf = uf.subs(MapPreSP);
        return uf;
    }
    
    void GPermutation(const ex & uf, const lst & xs) {
        auto expr = collect_ex(uf, xs);        
        map<ex,vector<ex>,ex_is_less> pgrp;
        
        if(true) { // only for polynomials
            auto cv_lst = collect_lst(expr, xs);
            exvector cvs;
            for(auto item : cv_lst) cvs.push_back(item);
            sort_vec(cvs); // first sort by coefficient
                    
            int nxi = xs.nops();
            bool first = true;
            lst xkey[nxi];
            lst subkey[nxi];
            ex clast;
            for(auto cv : cvs) {
                ex cc = cv.op(0);
                ex vv = cv.op(1);
                if(is_zero(cc)) continue;
                if(!first && !is_zero(cc-clast)) { // equal coefficient to the same key
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

            for(int i=0; i<nxi; i++) { // add last subkey
                sort_lst(subkey[i]); // need sort due to equal coefficient
                for(auto item : subkey[i]) xkey[i].append(item);
                subkey[i].remove_all();
            }
            
            exvector key_xi;
            for(int i=0; i<nxi; i++) key_xi.push_back(lst{xkey[i], xs.op(i)});
            sort_vec(key_xi); // first sort by key
            
            for(int i=0; i<nxi; i++) {
                auto xi = key_xi[i].op(1);
                auto xki = key_xi[i].op(0);
                pgrp[xki].push_back(xi);
            }
        }
 
        // pgrp - needs to permuation explicitly 
        long long npt = 1;
        for(auto pi : pgrp) {
            int nvi = pi.second.size();
            for(int i=1; i<=nvi; i++) npt *= i;  // nvi!
        }
        for(long long np=0; np<npt; np++) {
            long long npc = np;
            exmap xmap_tmp;
            for(auto pi : pgrp) {
                auto vi = pi.second;
                int nvi = vi.size();
                if(nvi<2) continue;
                long long npti = 1;
                for(int i=1; i<=nvi; i++) npti *= i; // nvi!
                int ck = npc % npti; // kth permutation
                npc = npc / npti;
            
                // https://stackoverflow.com/questions/1995328/are-there-any-better-methods-to-do-permutation-of-string
                auto vip = vi;
                int k = ck;
                for(int j=1; j<nvi; ++j) { // startoverflow: j starts from 1
                    std::swap(vip[k%(j+1)], vip[j]); 
                    k=k/(j+1);
                }
                
                for(int j=0; j<nvi; j++) if(vi[j]!=vip[j]) xmap_tmp[vi[j]] = vip[j];
            }
            if(xmap_tmp.size()>0) {
                ex expr_tmp = expr.subs(xmap_tmp,nopat);
                if(is_zero(expr-expr_tmp)) {
                    auto xs_tmp = xs.subs(xmap_tmp,nopat);
                    cout << xs_tmp << endl;
                }
            }
        }

    }

}
