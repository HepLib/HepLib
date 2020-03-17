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

    
    ex TIR(lst vis, lst eps) {
        lst lps;
        for(auto vi : vis) {
            if(!is_a<Pair>(vi) || !is_a<Vector>(vi.op(0)) || !is_a<Index>(vi.op(1)))
                throw Error("TIR invalid 1st arguments");
            lps.append(vi.op(0));
        }
        lps.sort();
        lps.unique();
        for(auto pi : eps) {
            if(!is_a<Vector>(pi)) throw Error("TIR invalid 2nd arguments");
        }
        
        auto visn = vis.nops();
        if(lps.nops()<2) {
            ex eqL=1, eqR=0;
            lst xs, is, bis, xs0;
            for(auto vi : vis) {
                eqL *= vi;
                is.append(vi.op(1));
            }
            
            for(int er=0; er<=vis.nops(); er=er+2) {
                auto ep_comb = combs(eps, er);
                for(auto item : ep_comb) {
                    ex bi=1;
                    for(int j=0; j<er; j++) bi *= SP(item.op(j), vis.op(j).op(1));
                    for(int j=er; j<visn; j=j+2) bi *= SP(vis.op(j).op(1), vis.op(j+1).op(1));
                    bi = bi.symmetrize(is);
                    bis.append(bi);
                    symbol x;
                    eqR += x*bi;
                    xs.append(x);
                    xs0.append(x==0);
            }}
            
            int n = bis.nops();
            lst eqns;
            for(auto bi : bis) {
                eqns.append((eqL+eqR).subs(sp_map)*bi);
            }
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
            for(int i=0; i<fvi; i++) ss << "(&J=fv" << i << ");" << endl;
            Fermat fermat;
            fermat.Init(InstallPrefix+"/bin/fer64");
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
            
            ss << "Array fM[" << n << "," << n+1 << "]" << endl;
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
cout << "mat2=" << mat2 << endl;
            ex res = 0;
            for(int i=0; i<n; i++) {
                res += bis.op(i) * mat2.op(i).op(n);
            }
            res = res.subs(sp_map);
            return res;
            
        }
        
        return 0;
        
    }

}

