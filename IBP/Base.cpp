/**
 * @file
 * @brief IBP
 */
 
#include "IBP.h"
#include <cmath>

namespace HepLib::IBP {

    namespace {
        ex_is_less less;
    }

    static void a_print(const ex & ex_in, const print_context & c) {
        c.s << "a[" << ex_in << "]";
    }

    REGISTER_FUNCTION(a, do_not_evalf_params().print_func<print_dflt>(a_print))
    
    void Base::Reduce() {
        Export();
        Run();
        Import();
    }
    
    void Base::Export(string garfn) {
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
    
    ex Base::TO() {
        lst shift;
        for(auto kv : Shift) shift.append(lst{kv.first, kv.second});
        
        return lst{ Internal, External, Replacements, Propagators, Cuts, DSP, ISP, shift, reCut, ProblemNumber, Symbol(WorkingDir), PIntegrals, MIntegrals, Rules, IsAlwaysZero };
    }
    
    void Base::Import(string garfn) {
        archive ar;
        ifstream in(garfn);
        in >> ar;
        in.close();
        
        Internal = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "Internal"));
        External = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "External"));
        Replacements = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "Replacements"));
        Propagators = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "Propagators"));
        Cuts = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "Cuts"));
        DSP = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "DSP"));
        ISP = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "ISP"));
        
        int n = ex_to<numeric>(ar.unarchive_ex(GiNaC_archive_Symbols, "NShift")).to_int();
        for(int i=0; i<n; i++) {
            int key = ex_to<numeric>(ar.unarchive_ex(GiNaC_archive_Symbols, ("ShiftK-"+to_string(i)).c_str())).to_int();
            ex val = ar.unarchive_ex(GiNaC_archive_Symbols, ("ShiftV-"+to_string(i)).c_str());
            Shift[key] = val;
        }
        
        reCut = !(ar.unarchive_ex(GiNaC_archive_Symbols,"reCut").is_zero()); // bool
        ProblemNumber = ex_to<numeric>(ar.unarchive_ex(GiNaC_archive_Symbols, "ProblemNumber")).to_int(); // int
        WorkingDir = ex_to<Symbol>(ar.unarchive_ex(GiNaC_archive_Symbols, "WorkingDir")).get_name();
        PIntegrals = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "PIntegrals"));
        MIntegrals = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "MIntegrals"));
        Rules = ex_to<lst>(ar.unarchive_ex(GiNaC_archive_Symbols, "Rules"));
        IsAlwaysZero = !(ar.unarchive_ex(GiNaC_archive_Symbols,"IsAlwaysZero").is_zero()); // bool
    
    }
    
    void Base::FROM(ex s) {
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
    
    /**
     * @brief Sort for all permuations, and return xs w.r.t. 1st permutation
     * @param in_expr the input expression, as the sort key, no need of polynormial of xs
     * @param xs the permutation list
     * @return 1st of sorted xs
     */
    lst SortPermutation(const ex & in_expr, const lst & xs) {
        auto expr = in_expr;
        bool isPoly = true;
        lst xRepl;
        map<ex,vector<int>,ex_is_less> pgrp;
        if(!expr.is_polynomial(xs)) {
            if(Verbose>2) cout << "SortPermutation: NOT a polynomials, try ALL permutations!" << endl;
            expr = expr.numer_denom();
            if(true || !expr.op(0).is_polynomial(xs) || !expr.op(1).is_polynomial(xs)) {
                isPoly = false;
                xRepl = xs;
                for(int i=0; i<xs.nops(); i++) pgrp[-1].push_back(i); // all permuations
            } else {
                expr = expr.op(0) * expr.op(1);
            }
        }
        
        if(isPoly) { // only for polynomials
            auto cv_lst = collect_lst(expr, xs);
            exvector cvs;
            for(auto item : cv_lst) cvs.push_back(item);
            sort_vec(cvs);
                    
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
            sort_vec(key_xi);
            
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
                ex expr2 = expr.subs(x2x);
                if(ex_less(expr2,expr1)) {
                    expr1 = expr2;
                    xRepl1 = ex_to<lst>(xRepl.subs(x2x));
                }
            }
            
            expr = expr1;
            xRepl = xRepl1;
        }
        return xRepl;
    }
    
    /**
     * @brief UF function
     * @param base the Base object
     * @param idx exponent for the internal Propagator
     * @return lst of {U, F, sign}
     */
    lst LoopUF(const Base & base, const ex & idx) {
        
        auto props = base.Propagators;
        
        // handle sign
        ex sign = 1;
        for(int i=0; i<props.nops(); i++) {
            auto ipr = props.op(i);
            
            if(ipr.has(iEpsilon)) {
                auto cc = ipr.coeff(iEpsilon);
                if(is_a<numeric>(cc)) {
                    if(cc<0) { // using +iEpsilon
                        sign = pow(-1, idx.op(i));
                        props.let_op(i) = ex(0)-props.op(i);
                    }
                    props.let_op(i) = props.op(i).subs(iEpsilon==0);
                    goto sign_done;
                } else throw Error("LoopUF: sign of iEpsilon NOT determined.");
            }
            
            for(auto lp : base.Internal) {
                if(ipr.degree(lp)==2) {
                    auto cc = ipr.coeff(lp,2);
                    if(is_a<numeric>(cc)) {
                        if(cc<0) { // using +l^2
                            sign = pow(-1, idx.op(i));
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
        exmap & cache = cache_by_prop[lst{props,base.Internal}];
        if(!using_cache || cache.find(key)==cache.end()) { // no cache item
            ut = 1;
            ft = expand_ex(ft);
            ft = subs_all(ft, base.Replacements);
            for(int i=0; i<base.Internal.nops(); i++) {
                auto t2 = ft.coeff(base.Internal.op(i),2);
                auto t1 = ft.coeff(base.Internal.op(i),1);
                auto t0 = ft.subs(base.Internal.op(i)==0);
                ut *= t2;
                if(is_zero(t2)) return lst{0,0,1};
                ft = exnormal(t0-t1*t1/(4*t2));
                ft = expand_ex(ft);
                ft = subs_all(ft, base.Replacements);
            }
            ft = exnormal(ut*ft);
            ft = exnormal(subs_all(ft, base.Replacements));
            ut = exnormal(subs_all(ut, base.Replacements));
            uf = exnormal(ut*ft);
            
            if(using_cache) cache[key] = lst{ut,ft,uf};
        } else {
            auto cc = cache[key];
            ut = cc.op(0);
            ft = cc.op(1);
            uf = cc.op(2);
        }

        ut = ut.subs(x2ax);
        ft = ft.subs(x2ax);
        uf = uf.subs(x2ax);
        
        uf = uf.subs(MapPreSP);
        auto xRepl = SortPermutation(uf,xs);
        for(int i=0; i<nxi; i++) xRepl.let_op(i)=(xRepl.op(i)==x(i));
        ut = (ut.subs(xRepl));
        ft = (ft.subs(xRepl));
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
            auto ipr = ps.op(i);
            
            if(ipr.has(iEpsilon)) {
                auto cc = ipr.coeff(iEpsilon);
                if(is_a<numeric>(cc)) {
                    if(cc<0) { // using +iEpsilon
                        sign = pow(-1, ns.op(i));
                        ps.let_op(i) = ex(0)-ps.op(i);
                    }
                    ps.let_op(i) = ps.op(i).subs(iEpsilon==0);
                    goto sign_done;
                } else throw Error("UF: sign of iEpsilon NOT determined.");
            }
            
            for(auto lp : loops) {
                if(ipr.degree(lp)==2) {
                    auto cc = ipr.coeff(lp,2);
                    if(is_a<numeric>(cc)) {
                        if(cc<0) { // using +l^2
                            sign = pow(-1, ns.op(i));
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
                            sign = pow(-1, ns.op(i));
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
                auto t0 = ft.subs(loops.op(i)==0);
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
                auto t0 = ft.subs(tloops.op(i)==0);
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
        ut1 = ut1.subs(x2ax);
        ut2 = ut2.subs(x2ax);
        ft = ft.subs(x2ax);
        uf = uf.subs(x2ax);
        
        uf = uf.subs(MapPreSP);
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
        return lst{ut1, ut2, ft, sign};
    }
    
    /**
     * @brief Find Rules for Integrals or Master Integrals
     * @param fs vector of Base pointer object
     * @param mi true for Master Integals
     * @param uf the function to compute the UF polynomial
     * @return rules replacement and left integrals or left master integrals
     */
    pair<exmap,lst> FindRules(vector<Base*> fs, bool mi, std::function<lst(const Base &, const ex &)> uf) {
        vector<pair<Base*,ex>> ibp_idx_vec;
        if(mi) for(auto fi : fs) for(auto item : fi->MIntegrals) ibp_idx_vec.push_back(make_pair(fi, item));
        else for(auto fi : fs) for(auto item : fi->Integrals) ibp_idx_vec.push_back(make_pair(fi, F(fi->ProblemNumber,item)));
        
        exvector uf_smi_vec = GiNaC_Parallel(ibp_idx_vec.size(), [&ibp_idx_vec,&uf](int idx)->ex {
            auto p = ibp_idx_vec[idx];
            const Base & fi = (*p.first);
            auto mi = p.second;
            auto ks = uf(fi,mi.subs(F(w1,w2)==w2));
            int nk = ks.nops()-1;
            lst key;
            for(int i=0; i<nk; i++) key.append(expand(ks.op(i)));
            lst pi = fi.Propagators; 
            sort_lst(pi);
            return lst{ key, lst{ pi, ks.op(nk), mi } }; 
            // ks.op(nk) -> the sign, 1st in lst used for sorting and must update n0=1 below
        }, "FR");
            
        map<ex,lst,ex_is_less> group;
        int ntotal = 0;
        for(auto item : uf_smi_vec) {
            ex key = item.op(0);
            group[key].append(item.op(1));
            ntotal++;
        }

        exmap rules;
        lst int_lst;
        for(auto g : group) {
            lst gs = ex_to<lst>(g.second);
            sort_lst(gs); // using sort_lst
            int n0 = 1; // update HERE
            auto c0 = gs.op(0).op(0+n0);
            auto v0 = gs.op(0).op(1+n0);
            for(int i=1; i<gs.nops(); i++) {
                auto ci = gs.op(i).op(0+n0);
                auto vi = gs.op(i).op(1+n0);
                rules[vi] = v0 * c0 / ci;
            }
            int_lst.append(v0);
        }
        
        if(Verbose>2) cout << "  \\--FindRules: " << WHITE << ntotal << " :> " << int_lst.nops() << RESET << " @ " << now(false) << endl;
        return make_pair(rules,int_lst);
    }
    
    static matrix Redrowech(matrix mat) {
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
                ss << mat(r,c).subs(iEpsilon==0).subs(v2f) << ",";
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
    
    bool Base::IsZero(ex sector) {
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
            auto t0 = ft.subs(Internal.op(i)==0);
            ut *= t2;
            if(is_zero(t2)) return true;
            ft = exnormal(t0-t1*t1/(4*t2));
            ft = expand(ft);
            ft = subs_all(ft, Replacements);
        }
        ut = exnormal(subs_all(ut, Replacements));
        ft = exnormal(ut*ft);
        ft = exnormal(subs_all(ft, Replacements));
    
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
            mat(ri,cn) = cv.op(0).subs(ks20);
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

}
