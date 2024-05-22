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
    
    ex UnContract(const ex expr, const lst &loop_ps, const lst &ext_ps) {
        // handle Eps/DGamma/Pair and related power
        int lproj = 0;
        return MapFunction([&lproj,loop_ps,ext_ps](const ex &e, MapFunction &self)->ex {
            string prefix = "dmi";
            int mode = 0;
            if(ext_ps.nops()>0) mode = 1;
            if(!Eps::has(e) && !DGamma::has(e) && (mode==0 || !Pair::has(e))) return e;
            else if(mode==1 && is_a<Pair>(e) && has_any(e,loop_ps) && has_any(e,ext_ps)) {
                Index idx(prefix+to_string(++lproj));
                auto p = ex_to<Pair>(e);
                return SP(p.op(0), idx) * SP(p.op(1), idx);
            } else if(is_a<Eps>(e)) {
                auto pis0 = ex_to<Eps>(e).pis;
                ex pis[4];
                ex cc = 1;
                for(int i=0; i<4; i++) {
                    pis[i] = pis0[i];
                    if(is_equal_any(pis[i],loop_ps)) {
                        Index idx(prefix+to_string(++lproj));
                        cc *= SP(pis[i], idx);
                        pis[i] = idx;
                    } else if(has_any(pis[i],loop_ps)) throw Error("UnContract: Eps still has loops.");
                }
                return LC(pis[0], pis[1], pis[2], pis[3]) * cc;
            } else if(is_a<DGamma>(e)) {
                Index idx(prefix+to_string(++lproj));
                auto g = ex_to<DGamma>(e);
                if(!is_equal_any(g.pi,loop_ps)) throw Error("UnContract: g.pi is NOT a loop.");
                return DGamma(idx, g.rl) * SP(g.pi, idx);
            } else if (e.match(TR(w))) {
                auto ret = self(e.op(0));
                auto cvs = collect_lst(ret, loop_ps);
                ret = 0;
                for(auto const & cv : cvs) {
                    if(DGamma::has(cv.op(1))) throw Error("UnContract: Not working for TIR.");
                    ret += TR(cv.op(0)) * cv.op(1);
                }
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
                if(!e.op(1).info(info_flags::posint)) return e; // no need to handle negative powers
                ex ret = 1;
                int pn = ex_to<numeric>(e.op(1)).to_int();
                for(int i=0; i<pn; i++) {
                    ret *= self(e.op(0));
                }
                return ret;
            } else return e.map(self);
        })(expr);
    }
    
    /**
     * @brief Tensor Index Reduction, note that we only handle numerator
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
        
        Fermat &fermat = Fermat::get();
        int &v_max = fermat.vmax;
        static exmap cache_map;
        
        auto expr = UnContract(expr_in, loop_ps); // UnContract
        auto cvs = collect_lst(expr, [&loop_ps](const ex & e)->bool{
            if(!Index::hasv(e)) return false;
            for(const_preorder_iterator i = e.preorder_begin(); i != e.preorder_end(); ++i) {
                auto item = *i;
                if(is_a<Pair>(item) && is_a<Index>(item.op(1)) && loop_ps.has(item.op(0))) return true;
            }
            return false;
        });
        
        expr = 0;
        for(auto cv : cvs) {
            auto const & cc = cv.op(0);
            auto const & vv = cv.op(1);
            if(vv.is_equal(1)) {
                expr += cc * vv;
                continue;
            }
            
            ex map_key = lst{vv, ext_ps};
            if(using_cache) {
                auto itr = cache_map.find(map_key);
                if(itr!=cache_map.end()) {
                    expr += cc * itr->second;
                    continue;
                }
            }
            lst vis, lps;
            map<ex,int,ex_is_less> pc;
            if(is_a<mul>(vv)) {
                for(auto item : vv) {
                    vis.append(item);
                    lps.append(item.op(0));
                    pc[item.op(0)]++;
                }
            } else {
                vis.append(vv);
                lps.append(vv.op(0));
            }
            lps.sort();
            lps.unique();
            
            ex res;
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
                ss.clear();
                ss.str("");
                
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

                res = 0;
                for(int i=0; i<n; i++) {
                    if(is_zero(mat2.op(i).op(i)) && !is_zero(mat2.op(i).op(n))) {
                        cout << cv << endl;
                        cout << mat << " :> " << mat2 << endl;
                        throw Error("No Solution in TIR.");
                    }
                    res += bis.op(i) * mat2.op(i).op(n);
                }
                res = res.subs(SP_map);
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
                for(auto lpi : lps) if(!is_zero(lpi-lp0)) ext_ps2.append(lpi);
                res = TIR(vv, lst{ lp0 }, ext_ps2);
                res = TIR(res, loop_ps, ext_ps);
            }
            
            res = exnormal(res);
            if(using_cache) cache_map.insert({map_key,res});
            expr += cc * res;
        }
    
        return expr;        
    }

}

