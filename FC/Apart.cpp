#include "FC.h"
#include "cln/cln.h"

namespace HepLib::FC {

    unsigned ApartIR1_SERIAL::serial = GiNaC::function::register_new(function_options("ApartIR",1).do_not_evalf_params().overloaded(2));
    unsigned ApartIR2_SERIAL::serial = GiNaC::function::register_new(function_options("ApartIR",2).do_not_evalf_params().overloaded(2));

    namespace {
        inline bool isOK(const ex &expr_in, const lst &vars) {
            ex item = expr_in;
            if(is_a<power>(item)) {
                if(!is_a<numeric>(item.op(1)) || !ex_to<numeric>(item.op(1)).is_integer()) return false;
                item = item.op(0);
            }
            bool isOK = true;
            for(auto v : vars) {
                if(item.degree(v)>1) return false;
            }
            return isOK;
        }
        
    }
    
    //--------------------------------------------------
    // ApartIR2ex
    //--------------------------------------------------
    ex ApartIR2ex(const ex & expr_in) {
        ex ret = expr_in;
        ret = MapFunction([](const ex & e, MapFunction &self)->ex{
            if(e.match(ApartIR(1)) || e.match(ApartIR(1,w))) return 1;
            else if(e.match(ApartIR(w1, w2))) {
                ex vars = e.op(1);
                ex ret=1;
                matrix mat = ex_to<matrix>(e.op(0));
                for(int c=0; c<mat.cols(); c++) {
                    ex sum=0;
                    for(int r=0; r<mat.rows()-2; r++) sum += mat(r,c) * vars.op(r);
                    sum += mat(mat.rows()-2,c);
                    ret *= pow(sum, mat(mat.rows()-1,c));
                }
                return ret;
            } else if(!e.has((ApartIR(w))) && !e.has((ApartIR(w1,w2)))) return e;
            else return e.map(self);
        })(ret);
        return ret;
    }
    
    //--------------------------------------------------
    // each column: [c1,...,cn,c0,n] -> (c1 x1+...+cn xn+c0)^n
    // return sum of coefficient * ApartIR
    //--------------------------------------------------
    ex Apart(const matrix & mat) {
        int nrow = mat.rows()-2;
        int ncol = mat.cols();
        lst null_vec;
        
        //--------------------------------------------------
        // null vector
        //--------------------------------------------------
        // if start feramt process many times, the efficient is slow
        bool use_fermat = false; 
        
        static exmap null_cache;
        if(is_zero(null_cache[sub_matrix(mat,0,nrow,0,ncol)])){
            if(use_fermat) {
                lst rep_vs;
                ex tree = mat;
                for(const_preorder_iterator i = tree.preorder_begin(); i != tree.preorder_end(); ++i) {
                    auto e = (*i);
                    if(is_a<symbol>(e) || is_a<Pair>(e) || is_a<Eps>(e)) {
                        rep_vs.append(e);
                    }
                }
                rep_vs.sort();
                rep_vs.unique();
                
                exmap v2f;
                symtab st;
                int fvi = 0;
                for(auto vi : rep_vs) {
                    auto name = "v" + to_string(fvi);
                    v2f[vi] = get_symbol(name);
                    st[name] = vi;
                    fvi++;
                }
                
                stringstream ss;
                for(int i=0; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
                Fermat fermat;
                fermat.Init();
                fermat.Execute(ss.str());
                ss.clear();
                ss.str("");
                
                ss << "Array m[" << nrow << "," << ncol+1 << "];" << endl;
                fermat.Execute(ss.str());
                ss.clear();
                ss.str("");
                
                ss << "[m]:=[(";
                for(int c=0; c<ncol; c++) {
                    for(int r=0; r<nrow; r++) {
                        ss << mat(r,c).subs(v2f) << ",";
                    }
                }
                for(int r=0; r<nrow; r++) ss << "0,";
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
                fermat.Exit();

                // make sure last char is 0
                if(ostr[ostr.length()-1]!='0') throw Error("TIR: last char is NOT 0.");
                ostr = ostr.substr(0, ostr.length()-1);
                string_trim(ostr);
                
                ostr.erase(0, ostr.find(":=")+2);
                string_replace_all(ostr, "[", "{");
                string_replace_all(ostr, "]", "}");
                Parser fp(st);
                auto mat2 = fp.Read(ostr);
                
                exmap xs;
                for(int c=0; c<ncol; c++) xs[c] = iWF(c);
                for(int r=nrow-1; r>=0; r--) {
                    ex xadd = 0;
                    int pi=-1;
                    for(int c=0; c<ncol; c++) {
                        if(!is_zero(get_op(mat2,r,c)) && pi<0) pi = c;
                        else xadd -= get_op(mat2,r,c)*xs[c];
                    }
                    xs[pi] = xadd;
                }
                for(int c=0; c<ncol; c++) null_vec.append(xs[c]);
            } else {
                matrix v(ncol, 1);
                lst sRepl;
                for (int c=0; c<ncol; c++) {
                    symbol t;
                    v(c,0) = t;
                    sRepl.append(t==iWF(c));
                }
                // Solve M*V = 0
                matrix zero(nrow, 1);
                matrix s = ex_to<matrix>(sub_matrix(mat,0,nrow,0,ncol)).solve(v,zero);
                for(int r=0; r<ncol; r++) null_vec.append(s(r,0).subs(sRepl));
            }
            null_cache[sub_matrix(mat,0,nrow,0,ncol)] = null_vec;
        } else {
            null_vec = ex_to<lst>(null_cache[sub_matrix(mat,0,nrow,0,ncol)]);
        }
        
        //--------------------------------------------------
        // check null & return ApartIR
        //--------------------------------------------------
        bool is_null = true;
        for(int c=0; c<ncol; c++) {
            if(!is_zero(null_vec.op(c))) {
                is_null = false;
                break;
            }
        }
        if(is_null) return ApartIR(mat);
        
        //--------------------------------------------------
        // handle numerator
        //--------------------------------------------------
        int ni=-1;
        for(int c=0; c<ncol; c++) {
            if(mat(nrow+1,c)>0 && !is_zero(null_vec.op(c))) {
                ni = c;
                break;
            }
        }
        
        if(ni!=-1) {
            auto nvec = subs(null_vec, iWF(w)==w+1);
            if(is_zero(nvec.op(ni))) nvec = subs(null_vec, iWF(w)==w*w+3);
            if(is_zero(nvec.op(ni))) throw Error("Apart: iWF to int failed with "+ex2str(null_vec.op(ni)));
            ex sol = 0;
            for(int c=0; c<ncol; c++) {
                if(c==ni) continue;
                sol -= nvec.op(c) * (iWF(c)-mat(nrow,c)); // iWF(c) refer to c-th column, and minus the const term
            }
            sol = sol/nvec.op(ni) + mat(nrow,ni); // last one: const term
            sol = mma_collect(pow(sol.subs(iEpsilon==0), mat(nrow+1,ni)), iWF(w)); // expand the numerator with power

            if(!is_a<add>(sol)) sol = lst{ sol };
            ex res = 0;
            for(auto item : sol) {
                int nzero = 0;
                matrix mat2(nrow+2, ncol-1);
                for(int c=0; c<ncol; c++) {
                    if(c==ni) continue;
                    int c2 = (c<ni ? c : c-1);
                    for(int r=0; r<nrow+1; r++) {
                        mat2(r,c2) = mat(r,c);
                    }
                    auto expn = mat(nrow+1,c) + item.degree(iWF(c));
                    if(is_zero(expn)) nzero++;
                    mat2(nrow+1,c2) = expn;
                }
                if(nzero>0) {
                    if((ncol-1-nzero)==0) {
                        res += ApartIR(1) * item.subs(iWF(w)==1);
                    } else {
                        matrix mat3(nrow+2, ncol-1-nzero);
                        int cc=0;
                        for(int c=0; c<ncol-1; c++) {
                            if(is_zero(mat2(nrow+1,c))) continue;
                            for(int r=0; r<nrow+2; r++) {
                                mat3(r,cc) = mat2(r,c);
                            }
                            cc++;
                        }
                        res += Apart(mat3) * item.subs(iWF(w)==1);
                    }
                } else {
                    res += Apart(mat2) * item.subs(iWF(w)==1);
                }
            }
            return res;
        }
        
        
        //--------------------------------------------------
        // handle all denominators
        //--------------------------------------------------
        
        ex cres0 = 0;
        for(int c=0; c<ncol; c++) cres0 += mat(nrow,c)*null_vec.op(c);
        cres0 = cres0.subs(iEpsilon==0);
        auto cres = subs(cres0, iWF(w)==w+1);
        auto nvec = subs(null_vec, iWF(w)==w+1);
        if(is_zero(cres)) {
            cres = subs(cres0, iWF(w)==w*w+3);
            nvec = subs(null_vec, iWF(w)==w*w+3);
        }

        // handle const is NOT zero
        if(!is_zero(cres)) {
            ex res=0;
            for(int c=0; c<ncol; c++) {
                if(is_zero(nvec.op(c))) continue;
                if(is_zero(mat(nrow+1,c)+1)) {
                    matrix mat2(nrow+2,ncol-1);
                    int ccc = 0;
                    for(int cc=0; cc<ncol; cc++) {
                        if(cc==c) continue;
                        for(int r=0; r<nrow+2; r++) mat2(r,ccc)=mat(r,cc);
                        ccc++;
                    }
                    res += Apart(mat2) * nvec.op(c);
                } else {
                    matrix mat2(mat);
                    mat2(nrow+1,c) = mat2(nrow+1,c)+1;
                    res += Apart(mat2) * nvec.op(c);
                }
            }
            res = res/cres;
            return res;
        } else {
            int ni=-1;
            for(int c=0; c<ncol; c++) {
                if(!is_zero(nvec.op(c))) {
                    ni = c;
                    break;
                }
            }

            ex res=0;
            for(int c=0; c<ncol; c++) {
                if(is_zero(nvec.op(c)) || c==ni) continue;
                if(is_zero(mat(nrow+1,c)+1)) {
                    matrix mat2(nrow+2,ncol-1);
                    int ccc = 0;
                    for(int cc=0; cc<ncol; cc++) {
                        if(cc==c) continue;
                        for(int r=0; r<nrow+2; r++) mat2(r,ccc)=mat(r,cc);
                        ccc++;
                    }
                    int ni2 = ni>c ? ni-1 : ni;
                    mat2(nrow+1,ni2) = mat2(nrow+1,ni2)-1;
                    res -= Apart(mat2) * nvec.op(c);
                } else {
                    matrix mat2(mat);
                    mat2(nrow+1,c) = mat2(nrow+1,c)+1;
                    mat2(nrow+1,ni) = mat2(nrow+1,ni)-1;
                    res -= Apart(mat2) * nvec.op(c);
                }
            }
            res = res/nvec.op(ni);
            return res;
        }
    }

    ex Apart(const ex &expr_in, const lst &vars_in, exmap sign_map) {
        exmap map1, map2;
        lst vars;
        for(auto v : vars_in) {
            symbol s;
            map1[v]=s;
            map2[s]=v;
            vars.append(s);
        }
        ex pref = 1;
        ex expr = expr_in.subs(map1);
        if(!is_a<mul>(expr)) expr = lst{expr};
        
        lst plst;
        for(auto item : expr) {
            bool has_var=false;
            for(auto v : vars) {
                if(item.has(v)) {
                    has_var=true;
                    break;
                }
            }
            if(has_var) {
                if(!isOK(item,vars)) throw Error("Apart: item is not linear wrt vars.");
                plst.append(item);
            } else pref *= item;
        }
        plst.sort();
        
        if(plst.nops()==0) return pref * ApartIR(1,vars_in);
        
        int nrow=vars.nops(), ncol=plst.nops();
        lst vars0;
        for(auto v : vars) vars0.append(v==0);
        matrix mat(nrow+2, ncol);
        for(int c=0; c<ncol; c++) {
            ex tmp = plst.op(c);
            if(is_a<power>(tmp)) {
                mat(nrow+1,c) = tmp.op(1);
                tmp = tmp.op(0);
            } else mat(nrow+1,c) = 1;
            
            // consider sign
            for(auto v : vars) {
                if(tmp.has(iEpsilon) && !is_zero(sign_map[iEpsilon])) { // iEpsilon first
                    ex sign = sign_map[iEpsilon]/tmp.coeff(iEpsilon);
                    pref /= pow(sign, mat(nrow+1,c));
                    tmp *= sign;
                    break;
                } else {
                    if(is_zero(tmp.coeff(v)) || is_zero(sign_map[v.subs(map2)])) continue;
                    ex sign = sign_map[v.subs(map2)]/tmp.coeff(v);
                    pref /= pow(sign, mat(nrow+1,c));
                    tmp *= sign;
                    break;
                }
            }
            
            for(int r=0; r<nrow; r++) {
                mat(r,c) = tmp.coeff(vars.op(r));
            }
            mat(nrow,c) = tmp.subs(vars0);
        }

        ex ret = Apart(mat);
        ret = MapFunction([vars](const ex & e, MapFunction &self)->ex{
            if(e.match(ApartIR(w))) return ApartIR(e.op(0), vars);
            else if(!e.has((ApartIR(w)))) return e;
            else return e.map(self);
        })(ret);
        ret = pref * ret;
        return ret.subs(map2);
    }
    
    ex Apart(const ex &expr_in, const lst &loops, const lst & extps) {
        auto expr = expr_in;
        lst sps;
        exmap sign;
        for(auto li : loops) {
            sign[SP(li)] = -1;
            for(auto li2: loops) {
                auto item = SP(li, li2).subs(sp_map);
                if(is_a<Pair>(item)) sps.append(item);
            }
            for(auto ei: extps) {
                auto item = SP(li, ei).subs(sp_map);
                if(is_a<Pair>(item)) sps.append(item);
            }
        }
        sps.sort();
        sps.unique();
        
        // TODO: maybe, to improve the numerator case
        if(expr.has(coVF(w))) throw Error("Apart: coVF exist.");
        expr = mma_collect(expr, loops, false, true);

        ex res = MapFunction([sps, sign](const ex & e, MapFunction &self)->ex{
            if(!e.has(coVF(w))) return e;
            else if(e.match(coVF(w))) return Apart(e.op(0), sps, sign);
            else return e.map(self);
        })(expr);

        // random check
        lst nlst;
        for(const_preorder_iterator i = res.preorder_begin(); i != res.preorder_end(); ++i) {
            auto e = (*i);
            if(is_a<symbol>(e) || is_a<Pair>(e)) nlst.append(e);
        }
        nlst.sort();
        nlst.unique();
        exmap nrepl;
        auto pi = cln::nextprobprime(3);
        for(auto ni : nlst) {
            pi = cln::nextprobprime(pi+1);
            nrepl[ni] = ex(1)/numeric(pi);
        }
        nrepl[iEpsilon]=0;
        ex chk = ApartIR2ex(subs(res,nrepl))-subs(expr_in,nrepl);
        chk = normal(chk);
        if(!is_zero(chk)) throw Error("Apart random check Failed.");
        
        return res;
    }
    
    ex ApartIRC(const ex & expr_in) {
        return MapFunction([](const ex & e, MapFunction &self)->ex {
            if(!e.has(ApartIR(w1,w2))) return e;
            else if(e.match(ApartIR(w1,w2))) {
                int n = e.op(1).nops();
                auto mat0 = ex_to<matrix>(e.op(0));
                int cc = mat0.cols();
                if(cc==n) return e;
                matrix mat(n+2,n);
                // zero each element
                for(int r=0; r<n+2; r++) {
                    for(int c=0; c<cc; c++) mat(r,c) = 0;
                }
                // n-row from mat0 to mat, note the last 2 rows still 0
                for(int r=0; r<n; r++) {
                    for(int c=0; c<cc; c++) mat(r,c) = mat0(r,c);
                }
                for(int i=0; i<n; i++) {
                    mat(i,cc) = 1;
                    auto r = mat.rank();
                    if(r==n) break;
                    if(r==cc+1) cc++;
                    else mat(i,cc) = 0;
                }
                if(mat.rank()!=n) throw Error("ApartIRC failed, NOT full rank.");
                // last 2 rows from mat0 to mat
                for(int r=n; r<n+2; r++) {
                    for(int c=0; c<mat0.cols(); c++) mat(r,c) = mat0(r,c);
                }
                return ApartIR(mat, e.op(1));
            } else return e.map(self);
        })(expr_in);
    }
    
    void Apart2FIRE(exvector &air_vec, lst vloops, lst vexts) {
        string wdir = to_string(getpid()) + "_FIRE";

        auto air_intg = 
        GiNaC_Parallel(-1, air_vec.size(), [air_vec] (int idx) {
            auto air = air_vec[idx];
            // fix ApartIR from archive
            air = air.subs(lst { 
                GiNaC::function(ApartIR1_SERIAL::serial, w1, w2)==ApartIR(w1,w2),
                GiNaC::function(ApartIR2_SERIAL::serial, w1)==ApartIR(w1)
            });
            
            air = ApartIRC(air);
            exset intg;
            find(air, ApartIR(w1, w2), intg);
            lst intgs;
            for(auto item : intg) intgs.append(item);
            return lst{air, intgs};
        }, "IRC");
        
        exset intg;
        for(int i=0; i<air_vec.size(); i++) {
            air_vec[i] = air_intg[i].op(0);
            for(auto item : air_intg[i].op(1)) intg.insert(item);
        }
        
        if(Verbose>0) cout << "  \\--Total Integrals: " << intg.size() << " @ " << now() << endl;
        
        exmap sp2;
        lst loops;
        for(auto vp1 : vloops) {
            auto p1 = ex_to<Vector>(vp1).name;
            loops.append(p1);
            for(auto vp2 : vloops) {
                auto p2 = ex_to<Vector>(vp2).name;
                sp2[SP(vp1,vp2)] = p1 * p2;
            }
            for(auto vp2 : vexts) {
                auto p2 = ex_to<Vector>(vp2).name;
                sp2[SP(vp1,vp2)] = p1 * p2;
            }
        }
        lst repls, exts;
        for(auto vp1 : vexts) {
            auto p1 = ex_to<Vector>(vp1).name;
            exts.append(p1);
            for(auto vp2 : vexts) {
                auto p2 = ex_to<Vector>(vp2).name;
                repls.append(p1*p2==SP(vp1,vp2).subs(sp_map));
            }
        }
        repls.sort();
        repls.unique();

        exmap ir2F;
        std::map<ex, FIRE*, ex_is_less> p2f;
        vector<FIRE*> fvec;
        int pn=1;
        for(auto item : intg) {
            auto mat = ex_to<matrix>(item.op(0));
            auto vars = ex_to<lst>(item.op(1));
            lst props, ns;
            int nrow = mat.rows();
            for(int c=0; c<mat.cols(); c++) {
                ex pc = 0;
                for(int r=0; r<nrow-2; r++) pc += mat(r,c) * vars.op(r);
                pc += mat(nrow-2,c);
                pc = pc.subs(sp2);
                if(Pair::has(pc)) {
                    cout << "sp2 = " << sp2 << endl;
                    throw Error("Apart2FIRE: Pair still exist: "+ex2str(pc));
                }
                props.append(pc);
                ns.append(ex(0)-mat(nrow-1,c).subs(sp2));
            }

            if(p2f[props]==NULL) {
                FIRE * f = new FIRE();
                p2f[props] = f;
                f->Propagators = props;
                f->Internal = loops;
                f->External = exts;
                f->Replacements = repls;
                f->WorkingDir = wdir;
                f->ProblemNumber = pn++;
                fvec.push_back(f);
            }
            FIRE * f = p2f[props];
            f->Integrals.append(ns);
            ir2F[item] = F(f->ProblemNumber, ns);
        }
        
        if(Verbose>0) cout << "  \\--Total FIRE: " << fvec.size() << " @ " << now() << endl;

        auto int_rules = FIRE::FindRules(fvec, false);

        for(auto &fp : fvec) {
            auto ints = fp->Integrals;
            for(int i=0; i<ints.nops(); i++) ints.let_op(i) = F(fp->ProblemNumber, ints.op(i));
            ints = ex_to<lst>(subs(ints, int_rules));
            ints.sort();
            ints.unique();
            for(int i=0; i<ints.nops(); i++) ints.let_op(i) = ints.op(i).op(1);
            fp->Integrals = ints;
        }

        auto fres= GiNaC_Parallel(-1, fvec.size(), [fvec](int idx)->ex {
            auto item = fvec[idx];
            item->Reduce();
            return lst {
                item->Dimension,
                item->MasterIntegrals,
                item->Rules
            };
        }, "FIRE");
        
        exmap F2F;
        for(int i=0; i<fres.size(); i++) {
            fvec[i]->Dimension = ex_to<numeric>(fres[i].op(0)).to_int();
            fvec[i]->MasterIntegrals = ex_to<lst>(fres[i].op(1));
            fvec[i]->Rules = ex_to<lst>(fres[i].op(2));
            for(auto item : fvec[i]->Rules) F2F[item.op(0)] = item.op(1);
        }
        
        MapFunction mapF([fvec](const ex &e, MapFunction &self)->ex {
            if(!e.has(F(w1,w2))) return e;
            else if(e.match(F(w1,w2))) {
                int idx = ex_to<numeric>(e.op(0)).to_int();
                return F(fvec[idx-1]->Propagators, e.op(1));
            } else return e.map(self);
        });
        
        auto mi_rules = FIRE::FindRules(fvec);

        for(auto &air : air_vec) {
            air = air.subs(ir2F);
            air = air.subs(int_rules);
            air = air.subs(F2F);
            air = air.subs(mi_rules);
            air = mapF(air);
        }
        
        for(auto fp : fvec) delete fp;
        system(("rm -rf "+wdir).c_str());
        
    }

}

