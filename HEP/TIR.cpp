/**
 * @file
 * @brief Functions to perform Tensor Index Reduction
 */
 
#include "HEP.h"
#include <cmath>

namespace HepLib {

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
            //for(int i=0; i<=n; i++) {
            for(int i=n; i>=0; i--) {
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
    
    /**
     * @brief Tensor Index Reduction
     * @param expr_in expression 
     * @param loop_ps lst contains loop vectors
     * @param ext_ps lst constains external vectors
     * @return TIR result
     */
    ex TIR(const ex &expr_in, const lst &loop_ps, const lst &ext_ps) {
        for(auto pi : loop_ps) {
            if(!is_a<Vector>(pi)) throw Error("TIR invalid 2nd argument");
        }
        for(auto pi : ext_ps) {
            if(!is_a<Vector>(pi)) throw Error("TIR invalid 3rd argument");
        }

        if(expr_in.has(coVF(w))) throw Error("TIR error: expr_in has coVF already.");
        
        // handle Eps/DGamma
        auto expr = expr_in;
        int lproj = 0;
        expr = MapFunction([&lproj,loop_ps](const ex &e, MapFunction &self)->ex {
            string prefix = "ti";
            if(!Eps::has(e) && !DGamma::has(e)) return e;
            else if(is_a<Eps>(e)) {
                auto pis0 = ex_to<Eps>(e).pis;
                ex pis[4];
                ex cc = 1;
                for(int i=0; i<4; i++) {
                    pis[i] = pis0[i];
                    if(is_equal_any(pis[i],loop_ps)) {
                        Index idx(prefix+to_string(++lproj));
                        cc *= SP(pis[i], idx);
                        pis[i] = idx;
                    } else if(has_any(pis[i],loop_ps)) throw Error("TIR: Eps still has loops.");
                }
                return LC(pis[0], pis[1], pis[2], pis[3]) * cc;
            } else if(is_a<DGamma>(e)) {
                Index idx(prefix+to_string(++lproj));
                auto g = ex_to<DGamma>(e);
                if(!is_equal_any(g.pi,loop_ps)) throw Error("TIR: g.pi is NOT a loop.");
                return DGamma(idx, g.rl) * SP(g.pi, idx);
            } else if (e.match(TR(w))) {
                auto ret = self(e.op(0));
                ret = collect_ex(ret, loop_ps, true);
                ret = ret.subs(coCF(w)==TR(w));
                return ret;
            } else if(is_a<add>(e)) {
                int lpj = lproj;
                int lpj_max = -100;
                ex ret = 0;
                for(auto item : e) {
                    lproj = lpj;
                    ret += self(item);
                    if(lpj_max<lproj) lpj_max = lproj;
                }
                lproj = lpj_max;
                return ret;
            } else if(is_a<power>(e)) {
                if(!e.op(1).info(info_flags::posint)) {
                    cout << e << endl;
                    throw Error("TIR: power is not info_flags::posint.");
                }
                ex ret = 1;
                int pn = ex_to<numeric>(e.op(1)).to_int();
                for(int i=0; i<pn; i++) {
                    ret *= self(e.op(0));
                }
                return ret;
            } else return e.map(self);
            throw Error("TIR: something should be wrong here");
            return 0;
        })(expr);

        expr = collect_ex(expr, [&loop_ps](const ex & e)->bool{
            if(!Index::hasv(e)) return false;
            for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
                auto item = *i;
                if(is_a<Pair>(item) && is_a<Index>(item.op(1)) && loop_ps.has(item.op(0))) return true;
            }
            return false;
        }, false, true);

        static exmap cache_map;
        expr = MapFunction([ext_ps,loop_ps](const ex &e, MapFunction &self)->ex{
            if(e.is_equal(coVF(1))) return 1;
            else if(!e.has(coVF(w))) return e;
            else if(e.match(coVF(w))) {
                ex map_key = lst{e,ext_ps};
                if(using_cache && cache_map.find(map_key)!=cache_map.end()) return cache_map[map_key];
                lst vis, lps;
                map<ex,int,ex_is_less> pc;
                if(is_a<mul>(e.op(0))) {
                    for(auto item : e.op(0)) {
                        vis.append(item);
                        lps.append(item.op(0));
                        pc[item.op(0)]++;
                    }
                } else {
                    vis.append(e.op(0));
                    lps.append(e.op(0).op(0));
                }
                lps.sort();
                lps.unique();
                
                auto visn = vis.nops();
                if(lps.nops()<2) {
                    ex eqL=1;
                    lst is, bis;
                    for(auto vi : vis) {
                        eqL *= vi;
                        is.append(vi.op(1));
                    }

                    //for(int er=((visn%2==0)?0:1); er<=visn; er+=2) {
                    for(int er=visn; er>=((visn%2==0)?0:1); er-=2) {
                        auto ep_comb = combs(ext_ps, er);
                        for(auto item : ep_comb) {
                            ex bi=1;
                            for(int j=0; j<er; j++) bi *= SP(item.op(j), is.op(j));
                            for(int j=er; j<visn; j=j+2) bi *= SP(is.op(j), is.op(j+1));
                            bi = bi.symmetrize(is);
                            bis.append(bi);
                    }}
                    
                    int n = bis.nops();
                    lst mat;
                    for(auto bi : bis) {
                        for(auto bj : bis) mat.append(bi*bj);
                    }
                    for(auto bj : bis) mat.append(eqL*bj);
                    
                    // we need to remap the dummy index, to avoid SP_map set the index object to 0
                    if(true) {
                        lst isu = is;
                        isu.sort();
                        isu.unique();
                        exmap i2u, u2i;
                        int c=0;
                        for(auto item : isu) {
                            auto ii = Index("TIR"+to_string(c++));
                            i2u[item] = ii;
                            u2i[ii] = item;
                        }
                        mat = ex_to<lst>(subs(mat,i2u));
                        mat = ex_to<lst>(form(mat).subs(u2i));
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
                    for(int i=0; i<fvi; i++) ss << "&(J=v" << i << ");" << endl;
                    Fermat fermat;
                    fermat.Init();
                    fermat.Execute(ss.str());
                    ss.clear();
                    ss.str("");
                    
                    ss << "Array m[" << n << "," << n+1 << "];" << endl;
                    fermat.Execute(ss.str());
                    ss.clear();
                    ss.str("");
                    
                    ss << "[m]:=[(";
                    for(auto item : mat) {
                        ss << item.subs(v2f) << ",";
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

                    ex res = 0;
                    for(int i=0; i<n; i++) {
                        if(is_zero(mat2.op(i).op(i)) && !is_zero(mat2.op(i).op(n))) throw Error("Zero Determinant in TIR.");
                        res += bis.op(i) * mat2.op(i).op(n);
                    }
                    res = res.subs(SP_map);
                    res = exnormal(res);
                    if(using_cache) cache_map[map_key] = res;
                    return res;
                } else {
                    int cmin=10000, cmax=-1;
                    ex pmin, pmax;
                    for(auto lpi : lps) {
                        auto cc = pc[lpi];
                        if(cc<cmin) { cmin = cc; pmin = lpi; }
                        if(cc>cmax) { cmax = cc; pmax = lpi; }
                    }
                    
                    auto lp0 = pmin;
                    auto ext_ps2 = ext_ps;
                    for(auto lpi : lps)
                        if(!is_zero(lpi-lp0)) ext_ps2.append(lpi);
                    ex ret = TIR(e.op(0), lst{ lp0 }, ext_ps2);
                    ret = TIR(ret, loop_ps, ext_ps);
                    ret = exnormal(ret);
                    if(using_cache) cache_map[map_key] = ret;
                    return ret;
                }
            } else return e.map(self);
        })(expr);
        
        return expr;        
    }

}

