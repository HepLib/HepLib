#include "FC.h"

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
        bool use_fermat = false;
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
            fermat.Init(InstallPrefix+"/bin/fer64");
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

            ostr = ostr.substr(0, ostr.length()-1);
            const char* WhiteSpace = " \t\v\r\n";
            if(!ostr.empty()) {
                ostr.erase(0, ostr.find_first_not_of(WhiteSpace));
                ostr.erase(ostr.find_last_not_of(WhiteSpace)+1);
            }
            
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
            if(is_zero(nvec.op(ni))) nvec = subs(null_vec, iWF(w)==w+3);
            if(is_zero(nvec.op(ni))) throw Error("Apart: iWF to int failed.");
            ex sol = 0;
            for(int c=0; c<ncol; c++) {
                if(c==ni) continue;
                sol -= nvec.op(c) * (iWF(c)-mat(nrow,c)); // iWF(c) refer to c-th column, and minus the const term
            }
            sol = sol/nvec.op(ni) + mat(nrow,ni); // last one: const term
            sol = mma_collect(pow(sol, mat(nrow+1,ni)), iWF(w)); // expand the numerator with power

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
        auto cres = subs(cres0, iWF(w)==w+1);
        auto nvec = subs(null_vec, iWF(w)==w+1);
        if(is_zero(cres)) {
            cres = subs(cres0, iWF(w)==w+3);
            nvec = subs(null_vec, iWF(w)==w+3);
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

    ex Apart(const ex &expr_in, const lst &vars) {
        ex pref = 1;
        ex expr = expr_in;
        if(!is_a<mul>(expr)) expr = lst{expr};
        
        lst plst;
        for(auto item : expr_in) {
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
        return pref * ret;
    }

}

