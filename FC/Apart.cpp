/**
 * @file
 * @brief Functions to perform partial fraction
 */

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
    
    /**
     * @brief convert ApartIR to ex
     * @param expr_in expression contains ApartIR
     * @return ApartIR converted into normal ex
     */
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
    
    /**
     * @brief convert ApartIR to F(ps, ns), ns is like FIRE convention
     * @param expr_in expression contains ApartIR
     * @return ApartIR converted into F(ps, ns)
     */
    ex ApartIR2F(const ex & expr_in) {
        ex ret = expr_in;
        ret = MapFunction([](const ex & e, MapFunction &self)->ex{
            if(e.match(ApartIR(1)) || e.match(ApartIR(1,w))) return 1;
            else if(e.match(ApartIR(w1, w2))) {
                ex vars = e.op(1);
                matrix mat = ex_to<matrix>(e.op(0));
                lst pns;
                for(int c=0; c<mat.cols(); c++) {
                    ex sum=0;
                    for(int r=0; r<mat.rows()-2; r++) sum += mat(r,c) * vars.op(r);
                    sum += mat(mat.rows()-2,c);
                    if(is_zero(mat(mat.rows()-1,c))) continue;
                    pns.append(lst{ sum, ex(0)-mat(mat.rows()-1,c) });
                }
                sort_lst(pns);
                lst ps, ns;
                for(auto item : pns) {
                    ps.append(item.op(0));
                    ns.append(item.op(1));
                }
                return F(ps, ns);
            } else if(!e.has((ApartIR(w))) && !e.has((ApartIR(w1,w2)))) return e;
            else return e.map(self);
        })(ret);
        return ret;
    }
    
    /**
     * @brief convert F(ps, ns) to normal ex, ns is like FIRE convention
     * @param expr_in expression contains F
     * @return F(ps, ns) converted into normal expression
     */
     ex F2ex(const ex & expr_in) {
        ex ret = expr_in;
        ret = MapFunction([](const ex & e, MapFunction &self)->ex{
            if(e.match(F(w1, w2))) {
                auto ps = e.op(0);
                auto ns = e.op(1);
                ex res = 1;
                for(int i=0; i<ps.nops(); i++) res *= pow(ps.op(i), ex(0)-ns.op(i));
                return res;
            } else if(!e.has((F(w1,w2)))) return e;
            else return e.map(self);
        })(ret);
        return ret;
     }
     
     /**
     * @brief convert F(ps, ns) to normal ex, ns is like FIRE convention
     * @param expr_in expression contains F
     * @return F(ps, ns) converted into normal expression wrapped in F
     */
     ex F2F(const ex & expr_in) {
        ex ret = expr_in;
        ret = MapFunction([](const ex & e, MapFunction &self)->ex{
            if(e.match(F(w1, w2))) {
                auto ps = e.op(0);
                auto ns = e.op(1);
                ex res = 1;
                for(int i=0; i<ps.nops(); i++) res *= pow(ps.op(i), ex(0)-ns.op(i));
                return F(res);
            } else if(!e.has((F(w1,w2)))) return e;
            else return e.map(self);
        })(ret);
        return ret;
     }
    
    /**
     * @brief Apart on matrix
     * @param mat each column: [c1,...,cn,c0,n] -> (c1 x1+...+cn xn+c0)^n
     * @return sum of coefficient * ApartIR
     */
    ex Apart(const matrix & mat) {
        int nrow = mat.rows()-2;
        int ncol = mat.cols();
        lst null_vec;
        
        // null vector
        bool use_fermat = true;
        static exmap null_cache;
        if(null_cache.find(sub_matrix(mat,0,nrow,0,ncol))==null_cache.end()){
            if(use_fermat) {
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
                    if(is_a<symbol>(e) || is_a<Pair>(e) || is_a<Eps>(e)) {
                        rep_vs.append(e);
                    }
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
                ss.clear();
                ss.str("");
                //fermat.Exit();
                
                // note the order, before exfactor (fermat_normal will be called again here)
                ss << "&(U=0);" << endl; // disable ugly printing
                ss << "@([m]);" << endl;
                ss << "&_G;" << endl;
                fermat.Execute(ss.str());
                ss.clear();
                ss.str("");

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
        
        // check null & return ApartIR
        bool is_null = true;
        for(int c=0; c<ncol; c++) {
            if(!is_zero(null_vec.op(c))) {
                is_null = false;
                break;
            }
        }
        if(is_null) return ApartIR(mat);
        
        // handle numerator
        int ni=-1;
        for(int c=0; c<ncol; c++) {
            if(mat(nrow+1,c)>0 && !is_zero(null_vec.op(c))) {
                ni = c;
                break;
            }
        }
        
        if(ni!=-1) {
            auto nvec = subs(null_vec, iWF(w)==w+1);
            if(is_zero(nvec.op(ni))) nvec = subs(null_vec, iWF(w)==w*w+1);
            if(is_zero(nvec.op(ni))) nvec = subs(null_vec, iWF(w)==w*w*w+1);
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
        
        // handle all denominators
        ex cres0 = 0;
        for(int c=0; c<ncol; c++) cres0 += mat(nrow,c)*null_vec.op(c);
        cres0 = cres0.subs(iEpsilon==0);
        auto cres = subs(cres0, iWF(w)==w+1);
        auto nvec = subs(null_vec, iWF(w)==w+1);
        if(is_zero(cres)) {
            cres = subs(cres0, iWF(w)==w*w+1);
            nvec = subs(null_vec, iWF(w)==w*w+1);
        }
        if(is_zero(cres)) {
            cres = subs(cres0, iWF(w)==w*w*w+1);
            nvec = subs(null_vec, iWF(w)==w*w*w+1);
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

    /**
     * @brief Apart on ex
     * @param expr_ino normal expresion, product of [ linear w.r.t. vars ]^n
     * @param vars_in independent variables
     * @param sgnmap the sign map
     * @return sum of coefficient * ApartIR
     */
    ex Apart(const ex &expr_ino, const lst &vars_in, exmap sgnmap) {
        // Apart on rational terms
        if(!is_a<lst>(expr_ino)) {
            auto cv_lst = mma_collect_lst(expr_ino, vars_in);
            ex res = 0;
            for(auto item : cv_lst) res += item.op(0) * Apart(lst{item.op(1)}, vars_in, sgnmap);
            
            // fermat_normal
            cv_lst = mma_collect_lst(res, ApartIR(w1,w2));
            res = 0;
            for(auto cv : cv_lst) res += fermat_normal(cv.op(0)) * cv.op(1);

            // random check
            lst nlst;
            for(const_preorder_iterator i = res.preorder_begin(); i != res.preorder_end(); ++i) {
                auto e = (*i);
                if(is_a<symbol>(e) || is_a<Pair>(e)) nlst.append(e);
            }
            for(auto var : vars_in) nlst.append(var);
            nlst.sort();
            nlst.unique();
            exmap nrepl;
            auto pi = cln::nextprobprime(3);
            for(auto ni : nlst) {
                pi = cln::nextprobprime(pi+1);
                nrepl[ni] = ex(1)/numeric(pi);
            }
            nrepl[iEpsilon]=0;
            ex chk = ApartIR2ex(subs(res,nrepl))-subs(expr_ino,nrepl);
            chk = normal(chk);
            if(!is_zero(chk)) throw Error("Apart@1 random check Failed.");
            return res;
        }
        
        // Apart on monomial term
        if(expr_ino.nops()!=1) throw Error("Apart: wrong convention found!");
        ex expr_in = expr_ino.op(0);
    
        exmap map1, map2;
        lst vars;
        for(int i=0; i<vars_in.nops(); i++) {
            auto v = vars_in.op(i);
            Symbol s("_apX"+to_string(i));
            map1[v]=s;
            map2[s]=v;
            vars.append(s);
        }
        
        ex expr = expr_in.subs(map1);
        if(!is_a<mul>(expr)) expr = lst{expr};
        
        // check only, try fermat_normal if faild
        bool ok = true;
        for(auto item : expr) {
            bool has_var=false;
            for(auto v : vars) {
                if(item.has(v)) {
                    has_var=true;
                    break;
                }
            }
            if(has_var) {
                if(!isOK(item,vars)) {
                    ok = false;
                    break;
                }
            }
        }
        if(!ok) {
            expr = expr_in.subs(map1);
            expr = fermat_normal(expr,true); // need option factor=true
            if(!is_a<mul>(expr)) expr = lst{expr};
        }
        
        lst pnlst;
        ex pref = 1;
        for(auto item : expr) {
            bool has_var=false;
            for(auto v : vars) {
                if(item.has(v)) {
                    has_var=true;
                    break;
                }
            }
            if(has_var) {
                if(!isOK(item,vars)) {
                    cout << expr_in << endl;
                    cout << item << endl;
                    throw Error("Apart: item is not linear wrt vars.");
                }
                
                ex pc, nc;
                if(is_a<power>(item)) {
                    pc = item.op(0);
                    nc = item.op(1);
                } else {
                    pc = item;
                    nc = 1;
                }
                
                // consider sign
                bool has_sgn = false;
                for(auto v : vars) {
                    if(pc.has(iEpsilon) && key_exists(sgnmap,iEpsilon)) { // iEpsilon first
                        ex sign = sgnmap[iEpsilon]/pc.coeff(iEpsilon);
                        pref /= pow(sign, nc);
                        pc *= sign;
                        has_sgn = true;
                        break;
                    } else {
                        if(is_zero(pc.coeff(v)) || !key_exists(sgnmap,v.subs(map2))) continue;
                        ex sign = sgnmap[v.subs(map2)]/pc.coeff(v);
                        pref /= pow(sign, nc);
                        pc *= sign;
                        has_sgn = true;
                        break;
                    }
                }
                if(!has_sgn) {
                    for(auto v : vars) {
                        if(is_zero(pc.coeff(v))) continue;
                        ex sign = 1/pc.coeff(v);
                        pref /= pow(sign, nc);
                        pc *= sign;
                        break;
                    }
                }
                pnlst.append(lst{ pc, nc });
            } else pref *= item;
        }
        sort_lst(pnlst);
        
        if(pnlst.nops()==0) return pref * ApartIR(1,vars_in);
        
        int nrow=vars.nops(), ncol=pnlst.nops();
        lst vars0;
        for(auto v : vars) vars0.append(v==0);
        matrix mat(nrow+2, ncol);
        for(int c=0; c<ncol; c++) {
            ex pn = pnlst.op(c);
            mat(nrow+1,c) = pn.op(1);
            auto tmp = pn.op(0);
            for(int r=0; r<nrow; r++) {
                mat(r,c) = tmp.coeff(vars.op(r));
            }
            mat(nrow,c) = tmp.subs(vars0);
        }

        ex ret = Apart(mat);
        auto cv_lst = mma_collect_lst(ret,ApartIR(w));
        ret = 0;
        for(auto cv : cv_lst) {
            ret += cv.op(0) * ApartIR(cv.op(1).op(0), vars);
        }
        ret = pref * ret;
        return ret.subs(map2);
    }
    
    /**
     * @brief Apart on ex
     * @param expr_ino input expression
     * @param loops list of loop Vector
     * @param extps list of external Vector
     * @return sum of coefficient * ApartIR
     */
    ex Apart(const ex &expr_ino, const lst &loops, const lst & extps) {
        auto expr_in = expr_ino.subs(SP_map);
        auto expr = expr_in;
        
        lst sps;
        exmap sgnmap;
        for(auto li : loops) {
            sgnmap[SP(li)] = 1;
            for(auto li2: loops) {
                auto item = SP(li, li2).subs(SP_map);
                if(is_a<Pair>(item)) sps.append(item);
            }
            for(auto ei: extps) {
                auto item = SP(li, ei).subs(SP_map);
                if(is_a<Pair>(item)) sps.append(item);
            }
        }
        sps.sort();
        sps.unique();
        sort_lst(sps);
        
        auto cv_lst = mma_collect_lst(expr, loops);
        ex res = 0;
        for(auto item : cv_lst) {
            res += item.op(0) * Apart(lst{item.op(1)}, sps, sgnmap);
        }
        
        // fermat_normal
        cv_lst = mma_collect_lst(res, ApartIR(w1,w2));
        res = 0;
        for(auto cv : cv_lst) res += fermat_normal(cv.op(0)) * cv.op(1);

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
        if(!is_zero(chk)) throw Error("Apart@2 random check Failed.");
        
        return res;
    }
    
    /**
     * @brief complete the ApartIR elements
     * @param expr_in input expression
     * @return ApartIR with complete matrix rank, ready for IBP reduction
     */
    ex ApartIRC(const ex & expr_in) {
        return MapFunction([](const ex & e, MapFunction &self)->ex {
            if(!e.has(ApartIR(w1,w2))) return e;
            else if(e.match(ApartIR(w1,w2))) {
                int n = e.op(1).nops();
                auto mat0 = ex_to<matrix>(e.op(0));
                matrix mat(n+2,n);
                int cc = mat0.cols();
                if(cc==n) mat=mat0;
                else {
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
                }
                return ApartIR(mat, e.op(1));
            } else return e.map(self);
        })(expr_in);
    }
    
    /**
     * @brief perform IBP reduction on the Aparted input
     * @param IBPmethod ibp method used, 0-No IBP, 1-FIRE, 2-KIRA
     * @param air_vec vector contains aparted input, ApartIRC will be call internally 
     * @param loops_exts lst { loop vectors, external vectors, }
     * @param cut_props cut propagators, default is { }
     * @param uf the function to compute UF polynomial
     * @return nothing returned, the input air_vec will be updated
     */
    void ApartIBP(int IBPmethod, exvector &air_vec, const lst & loops_exts, const lst & cut_props, 
        std::function<lst(const Base &, const ex &)> uf) {
        
        if(loops_exts.nops()<2) throw Error("loops_exts size() < 2;");
        
        string wdir = to_string(getpid());
        if(IBPmethod==1) wdir = wdir + "_FIRE";
        else if(IBPmethod==2) wdir = wdir + "_KIRA";
        
        bool aparted = false;
        for(auto air : air_vec) {
            if(air.has(ApartIR(w1,w2))) {
                aparted = true;
                break;
            }
        }
        
        auto air_intg = 
        GiNaC_Parallel(air_vec.size(), [aparted,air_vec,loops_exts] (int idx) {
            auto air = air_vec[idx];
            if(!aparted) {
                lst lmom = ex_to<lst>(loops_exts.op(0));
                lst emom = ex_to<lst>(loops_exts.op(1));
                air = Apart(air,lmom,emom);
            }
            air = air.subs(SP_map);            
            air = ApartIRC(air);
            exset intg;
            find(air, ApartIR(w1, w2), intg);
            lst intgs;
            for(auto item : intg) intgs.append(item);
            return lst{air, intgs};
        }, "Apart");
        
        exvector intg_sort;
        if(true) {
            exset intg;
            for(int i=0; i<air_vec.size(); i++) {
                air_vec[i] = air_intg[i].op(0);
                for(auto item : air_intg[i].op(1)) intg.insert(item);
            }
            air_intg.clear();
            air_intg.shrink_to_fit();
            for(auto item : intg) intg_sort.push_back(item);
            sort_vec(intg_sort);
        }
                
        lst repls;
        auto sps = sp_map();
        for(auto kv : sps) repls.append(kv.first == kv.second);
        
        lst loops, exts;
        for(auto li : loops_exts.op(0)) {
            if(is_a<Vector>(li)) loops.append(ex_to<Vector>(li).name);
            else loops.append(li);
        }
        for(auto li : loops_exts.op(1)) {
            if(is_a<Vector>(li)) exts.append(ex_to<Vector>(li).name);
            else exts.append(li);
        }
        
        exmap IR2F;
        std::map<ex, IBP::Base*, ex_is_less> p2IBP;
        vector<IBP::Base*> ibp_vec;
        int pn=1;
        for(auto item : intg_sort) {
            auto mat = ex_to<matrix>(item.op(0));
            auto vars = ex_to<lst>(item.op(1));
            lst pns;
            int nrow = mat.rows();
            for(int c=0; c<mat.cols(); c++) {
                ex pc = 0;
                for(int r=0; r<nrow-2; r++) pc += mat(r,c) * vars.op(r);
                pc += mat(nrow-2,c);
                pc = SP2sp(pc);
                pns.append(lst{ pc, ex(0)-mat(nrow-1,c) }); // note the convension
            }
            sort_lst(pns);
            for(auto cut : cut_props) pns.prepend(lst{ SP2sp(item.subs(SP_map)), 1 });
            
            lst props, ns;
            for(auto item : pns) {
                props.append(item.op(0));
                ns.append(item.op(1));
            }

            if(p2IBP[props]==NULL) {
                IBP::Base* ibp;
                if(IBPmethod==0) ibp = new Base();
                else if(IBPmethod==1) ibp = new FIRE();
                else if(IBPmethod==2) ibp = new KIRA();
                else {
                    ibp = new Base();
                    IBPmethod = 0;
                }
                
                p2IBP[props] = ibp;
                ibp->Propagators = props;
                ibp->Internal = loops;
                ibp->External = exts;
                ibp->Replacements = repls;
                ibp->Pairs = ex_to<lst>(SP2sp(Pair::all(vars)));
                ibp->WorkingDir = wdir;
                ibp->ProblemNumber = pn++;
                if(cut_props.nops()>0) {
                    for(int i=0; i<cut_props.nops(); i++) ibp->Cuts.append(i+1);
                }
                ibp_vec.push_back(ibp);
            }
            IBP::Base* ibp = p2IBP[props];
            ibp->Integrals.append(ns);
            IR2F[item] = F(ibp->ProblemNumber, ns);
        }
        
        if(Verbose>0) cout << "  \\--Total Ints/Pros: " << intg_sort.size() << "/" << ibp_vec.size() << " @ " << now(false) << endl;
        
        vector<Base*> base_vec;
        for(auto ibp : ibp_vec) base_vec.push_back(ibp);
        auto rules_ints = FindRules(base_vec, false, uf);

        map<int,lst> pn_ints_map;
        for(auto item : rules_ints.second) {
            int pn = ex_to<numeric>(item.op(0)).to_int();
            pn_ints_map[pn].append(item.op(1));
        }

        vector<IBP::Base*> ibp_vec_re;
        int nints = 0;
        for(auto pi : pn_ints_map) {
            auto ibp = ibp_vec[pi.first-1];
            ibp->Integrals = pi.second;
            nints += ibp->Integrals.nops();
            ibp_vec_re.push_back(ibp);
        }

        if(Verbose>0) cout << "  \\--Refined Ints/Pros: " << nints << "/" << ibp_vec_re.size() << " @ " << now(false) << endl;
        
        MapFunction _F2ex([ibp_vec](const ex &e, MapFunction &self)->ex {
            if(!e.has(F(w1,w2))) return e;
            else if(e.match(F(w1,w2))) {
                int idx = ex_to<numeric>(e.op(0)).to_int();
                auto pso = ex_to<lst>(ibp_vec[idx-1]->Propagators);
                auto nso = ex_to<lst>(e.op(1));
                lst ps, ns;
                for(int i=0; i<pso.nops(); i++) {
                    if(nso.op(i).is_zero()) continue;
                    ps.append(pso.op(i));
                    ns.append(nso.op(i));
                }
                return F(ps,ns);
            } else return e.map(self);
        });
        
        if(IBPmethod==0) {
            auto air_res =
            GiNaC_Parallel(air_vec.size(), 1, [&](int idx)->ex {
                auto air = air_vec[idx];
                air = subs_naive(air,IR2F);
                air = subs_naive(air,rules_ints.first);
                air = _F2ex(air);
                auto cv_lst = mma_collect_lst(air, F(w1,w2));
                air = 0;
                for(auto cv : cv_lst) air += cv.op(1) * fermat_normal(cv.op(0));
                return air;
            }, "A2F");
            
            for(auto fp : ibp_vec) delete fp;
            system(("rm -rf "+wdir).c_str());

            for(int i=0; i<air_vec.size(); i++) air_vec[i] = air_res[i];   
            return;
        }
        
        if(IBPmethod==1) {
            for(auto ibp : ibp_vec_re) ibp->Export();
            auto nproc = 2*CpuCores()/FIRE::Threads;
            int cproc = 0;
            if(nproc<2) nproc = 2;
            #pragma omp parallel for num_threads(nproc) schedule(dynamic, 1)
            for(int pi=0; pi<ibp_vec_re.size(); pi++) {
                if(Verbose>1) {
                    #pragma omp critical
                    cout << "\r                                        \r" << "  \\--FIRE Reduction [" << (++cproc) << "/" << ibp_vec_re.size() << "] " << flush;
                }
                ibp_vec_re[pi]->Run();
            }
            if(Verbose>1) cout << "@" << now(false) << endl;
            for(auto item : ibp_vec_re) item->Import();
            system(("rm -rf "+wdir).c_str());
        } else if(IBPmethod==2) {
            for(auto ibp : ibp_vec_re) ibp->Reduce();
            system(("rm -rf "+wdir).c_str());
        }
        
        vector<Base*> base_re;
        for(auto f : ibp_vec_re) base_re.push_back(f);
        auto mi_rules = FindRules(base_re, true, uf);
                        
        auto rules_vec =
        GiNaC_Parallel(ibp_vec_re.size(), 10, [&](int idx)->ex {
            lst rules;
            for(auto item : ibp_vec_re[idx]->Rules) {
                auto rr = item.op(1);
                rr = subs_naive(rr,mi_rules.first);
                auto cv_lst = mma_collect_lst(rr, F(w1,w2));
                rr = 0;
                for(auto cv : cv_lst) rr += cv.op(1) * fermat_normal(cv.op(0));
                rules.append(item.op(0)==rr);
            }
            return rules;
        }, "F2F");
                
        exmap F2F;
        for(auto rule : rules_vec) {
            for(auto lr : rule) F2F[lr.op(0)] = lr.op(1);
        }
        
        auto air_res =
        GiNaC_Parallel(air_vec.size(), 1, [&](int idx)->ex {
            auto air = air_vec[idx];
            air = subs_naive(air,IR2F);
            air = subs_naive(air,rules_ints.first);
            air = subs_naive(air,mi_rules.first);
            air = subs_naive(air,F2F);
            air = _F2ex(air);
            auto cv_lst = mma_collect_lst(air, F(w1,w2));
            air = 0;
            for(auto cv : cv_lst) air += cv.op(1) * fermat_normal(cv.op(0));
            return air;
        }, "A2F");
            
        for(auto fp : ibp_vec) delete fp;

        for(int i=0; i<air_vec.size(); i++) air_vec[i] = air_res[i];        
    }

}
