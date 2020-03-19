#include "FC.h"
#include <cmath>

namespace HepLib::FC {

    namespace {
        
        lst combs(const lst & items, int n) {
            lst ret;
            if(n<1 || items.nops()<1) return lst{ret};
            if(n==1) {
                for(auto it : items) ret.append(lst{it});
                return ret;
            }
            if(items.nops()==1) {
                lst comb;
                for(int i=0; i<n; i++) comb.append(items.op(0));
                ret.append(comb);
                return ret;
            }
            auto its = items;
            ex first = its.op(0);
            its.remove_first();
            for(int i=0; i<=n; i++) {
                auto rets = combs(its, n-i);
                for(auto it : rets) {
                    lst item = ex_to<lst>(it);
                    for(int j=0; j<i; j++) item.append(first);
                    ret.append(item);
                }
            }
            return ret;
        }
        
    }

    ex TIR(const ex &expr_in, const lst &loop_ps, const lst &ext_ps) {
        for(auto pi : loop_ps) {
            if(!is_a<Vector>(pi)) throw Error("TIR invalid 2nd argument");
        }
        for(auto pi : ext_ps) {
            if(!is_a<Vector>(pi)) throw Error("TIR invalid 3rd argument");
        }
        
        if(expr_in.has(coVF(w))) throw Error("TIR error: expr_in has coVF already.");
        auto expr = mma_collect(expr_in, [&](const ex & e)->bool{
            for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
                auto item = *i;
                if(is_a<Pair>(item) && is_a<Index>(item.op(1)) && loop_ps.has(item.op(0))) return true;
            }
            return false;
        }, false, true);
        
        expr = MapFunction([&](const ex &e, MapFunction &self)->ex{
            if(e.match(coVF(w))) {
                lst vis, lps;
                if(is_a<mul>(e.op(0))) {
                    for(auto item : e.op(0)) {
                        vis.append(item);
                        lps.append(item.op(0));
                    }
                } else {
                    vis.append(e.op(0));
                    lps.append(e.op(0).op(0));
                }
                lps.sort();
                lps.unique();
                
                auto visn = vis.nops();
                if(lps.nops()<2) {
                    ex eqL=1, eqR=0;
                    lst xs, is, bis, xs0;
                    for(auto vi : vis) {
                        eqL *= vi;
                        is.append(vi.op(1));
                    }
                    
                    for(int er=((visn%2==0) ? 0 : 1); er<=visn; er=er+2) {
                        auto ep_comb = combs(ext_ps, er);
                        for(auto item : ep_comb) {
                            ex bi=1;
                            for(int j=0; j<er; j++) bi *= SP(item.op(j), vis.op(j).op(1));
                            for(int j=er; j<visn; j=j+2) bi *= SP(vis.op(j).op(1), vis.op(j+1).op(1));
                            bi = bi.symmetrize(is);
                            bi = bi.normal().numer(); // drop denominator
                            bis.append(bi);
                            symbol x;
                            eqR += x*bi;
                            xs.append(x);
                            xs0.append(x==0);
                    }}
                    
                    vector<ex> bis2;
                    for(auto item : bis) bis2.push_back(item);
                    sort(bis2.begin(), bis2.end(), [](const ex &a, const ex &b)->bool {
                        int an = 1, bn = 1;
                        if(is_a<add>(a)) an = a.nops();
                        if(is_a<add>(b)) bn = b.nops();
                        return an<bn;
                    });
                    bis.remove_all();
                    for(auto item : bis2) bis.append(item);
                    
                    int n = bis.nops();
                    lst eqns;
                    for(auto bi : bis) eqns.append((eqL+eqR).subs(sp_map)*bi);
                    eqns = ex_to<lst>(form(eqns).subs(sp_map));

                    matrix mat(n, n+1);
                    for(int i=0; i<n; i++) {
                        auto eqn = eqns.op(i);
                        for(int j=0; j<n; j++) {
                            mat(i,j) = eqn.coeff(xs.op(j));
                        }
                        mat(i, n) = eqn.subs(xs0);
                    }
                    
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
                        auto name = "fv" + to_string(fvi);
                        v2f[vi] = get_symbol(name);
                        st[name] = vi;
                        fvi++;
                    }

                    stringstream ss;
                    for(int i=0; i<fvi; i++) ss << "&(J=fv" << i << ");" << endl;
                    Fermat fermat;
                    fermat.Init(InstallPrefix+"/bin/fer64");
                    fermat.Execute(ss.str());
                    ss.clear();
                    ss.str("");
                    
                    ss << "Array fM[" << n << "," << n+1 << "];" << endl;
                    fermat.Execute(ss.str());
                    ss.clear();
                    ss.str("");

                    ss << "[fM]:=[(";
                    for(int c=0; c<n+1; c++) {
                        for(int r=0; r<n; r++) {
                            ss << mat(r,c).subs(v2f);
                            if(!(c==n && r==n-1)) ss << ",";
                        }
                    }
                    ss << ")];" << endl;
                    ss << "Redrowech([fM]);" << endl;
                    fermat.Execute(ss.str());
                    ss.clear();
                    ss.str("");
                    
                    ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
                    ss << "![fM" << endl;
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

                    ex res = 0;
                    for(int i=0; i<n; i++) {
                        if(is_zero(mat2.op(i).op(i))) throw Error("Zero Determinant in TIR.");
                        res += bis.op(i) * mat2.op(i).op(n);
                    }
                    res = res.subs(sp_map);
                    return res;
                } else {
                    auto lp1st = lps.op(0);
                    auto ext_ps2 = ext_ps;
                    for(int i=1; i<lps.nops(); i++) ext_ps2.append(lps.op(i));
                    ex ret = TIR(e.op(0), lst{lp1st}, ext_ps2);
                    return TIR(ret, loop_ps, ext_ps2);
                }
            } else return e.map(self);
        })(expr);
    
        return expr;
        
    }

}

