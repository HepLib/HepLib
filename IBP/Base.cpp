/**
 * @file
 * @brief IBP
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib::IBP {

    static void a_print(const ex & ex_in, const print_context & c) {
        c.s << "a[" << ex_in << "]";
    }

    REGISTER_FUNCTION(a, do_not_evalf_params().print_func<print_dflt>(a_print))
    
    void Base::Reduce() {
        Export();
        Run();
        Import();
    }
    
    /**
     * @brief Sort for all permuations, and return xs w.r.t. 1st permutation
     * @param expr the input expression, as the sort key, no need of polynormial of xs
     * @param xs the permutation list
     * @return 1st of sorted xs
     */
    lst SortPermutation(const ex & in_expr, const lst & xs) {
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
            auto cv_lst = mma_collect_lst(expr, xs);
            exvector cvs;
            for(auto item : cv_lst) cvs.push_back(item);
            sort_vec_by(cvs,0);
                    
            int nxi = xs.nops();
            bool first = true;
            lst xkey[nxi];
            lst subkey[nxi];
            ex clast;
            for(auto cv : cvs) {
                ex cc = cv.op(0);
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
            
            // https://stackoverflow.com/questions/1995328/are-there-any-better-methods-to-do-permutation-of-string
            long long nf = 1;
            for(int i=2; i<=nvi; i++) nf *= i;
            for(long long ck=1; ck<=nf; ck++) { // n! permutations
                auto vi2 = vi;
                int k = ck;
                for(int j=1; j<nvi; ++j) {
                    std::swap(vi2[k%(j+1)], vi2[j]); 
                    k=k/(j+1);
                }
                
                exmap x2x;
                for(int j=0; j<nvi; j++) {
                    if(vi[j]!=vi2[j]) x2x[xs.op(vi[j])]=xs.op(vi2[j]);
                }
                ex expr2 = subs_naive(expr,x2x);
                if(ex_less(expr2,expr1)) {
                    expr1 = expr2;
                    xRepl1 = ex_to<lst>(subs_naive(xRepl,x2x));
                }
            }
            
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
    lst LoopUF(const Base & fire, const ex & idx) {
        static map<ex,exmap,ex_is_less> cache_by_prop;
        exmap & cache = cache_by_prop[lst{fire.Propagators,fire.Internal}];
        ex ut, ft, uf;
        lst key;
        lst xs;
        exmap x2ax;
        int nxi=0;
        int nps = fire.Propagators.nops();
        ft = 0;
        for(int i=0; i<nps; i++) {
            if(is_zero(idx.op(i))) {
                key.append(0);
                continue;
            }
            key.append(1);
            if(!is_zero(idx.op(i)-1)) x2ax[x(nxi)] = a(idx.op(i)) * x(nxi);
            xs.append(x(nxi));
            ft -= x(nxi) * fire.Propagators.op(i); // only used when no cache item
            nxi++;
        }
        
        if(cache.find(key)==cache.end()) { // no cache item
            ut = 1;
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
            ft = normal(ut*ft);
            ft = normal(subs_all(ft, fire.Replacements));
            ut = normal(subs_all(ut, fire.Replacements));
            uf = normal(ut*ft);
            
            cache[key] = lst{ut,ft,uf};
        } else {
            auto cc = cache[key];
            ut = cc.op(0);
            ft = cc.op(1);
            uf = cc.op(2);
        }
        ut = ut.subs(x2ax);
        ft = ft.subs(x2ax);
        uf = uf.subs(x2ax);
        
        auto xRepl = SortPermutation(uf,xs);
        for(int i=0; i<nxi; i++) xRepl.let_op(i)=(xRepl.op(i)==x(i));
        ut = (subs_naive(ut,xRepl)); 
        ft = (subs_naive(ft,xRepl));
        return lst{ut, ft};
    }  
    
    /**
     * @brief UF function, from FIRE.m
     * @param idx exponent for the internal Propagator
     * @return lst of {U, F}
     */
    lst UF(const ex & ps, const ex & ns, const ex & loops, const ex & tloops, const ex & lsubs, const ex & tsubs) {
        static map<ex,exmap,ex_is_less> cache_by_prop;
        exmap & cache = cache_by_prop[lst{ps,loops,tloops}];
        
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
        
        if(cache.find(key)==cache.end()) { // no cache item
            ut1 = 1;
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
            ft = normal(ut1*ft);
            ft = normal(subs_all(ft, lsubs));
            ut1 = normal(subs_all(ut1, lsubs));

            ut2 = 1;
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
            ft = normal(ut2*ft);
            ft = normal(subs_all(ft, tsubs));
            ut2 = normal(subs_all(ut2, tsubs));
            
            uf = normal(ut1*ut2*ft);
            cache[key] = lst{ut1,ut2,ft,uf};
        } else {
            auto cc = cache[key];
            ut1 = cc.op(0);
            ut2 = cc.op(1);
            ft = cc.op(2);
            uf = cc.op(3);
        }
        ut1 = ut1.subs(x2ax);
        ut2 = ut2.subs(x2ax);
        ft = ft.subs(x2ax);
        uf = uf.subs(x2ax);
        
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

        ut1 = (ut1.subs(xRepl));
        ut2 = (ut2.subs(xRepl));
        ft = (ft.subs(xRepl));
        return lst{ut1, ut2, ft};
    }
    
    /**
     * @brief Find Rules for Integrals or Master Integrals
     * @param fs vector of Base pointer object
     * @param mi true for Master Integals
     * @return rules replacement and left integrals or left master integrals
     */
    pair<exmap,lst> FindRules(vector<Base*> & fs, bool mi, std::function<lst(const Base &, const ex &)> uf) {
        exvector uf_mi_vec;
        if(mi) {
            uf_mi_vec = GiNaC_Parallel(fs.size(), [mi,fs,uf](int idx)->ex {
                const Base & fi = *(fs[idx]); // only here
                lst uf_mi_lst;
                for(auto mi : fi.MasterIntegrals) {
                    uf_mi_lst.append(lst{ uf(fi,mi.subs(F(w1,w2)==w2)), mi });
                }
                return uf_mi_lst;
            }, "MRules");
        } else {
            uf_mi_vec = GiNaC_Parallel(fs.size(), [mi,fs,uf](int idx)->ex {
                const Base & fi = *(fs[idx]); // only here
                lst uf_mi_lst;
                for(auto mi : fi.Integrals) {
                    uf_mi_lst.append(lst{ uf(fi,mi), F(fi.ProblemNumber,mi) });
                }
                return uf_mi_lst;
            }, "FRules");
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
            sort_lst(gs);
            for(int i=1; i<gs.nops(); i++) rules[gs.op(i)]=gs.op(0);
            int_lst.append(gs.op(0));
        }
        
        if(Verbose>2) cout << "  \\--FindRules: " << ntotal << " :> " << int_lst.nops() << " @ " << now(false) << endl;
        return make_pair(rules,int_lst);
    }

}
