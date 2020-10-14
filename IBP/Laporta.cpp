#include "Process.h"
#include "IBP.h"

namespace HepLib::IBP {

lst Laporta::F2lst(ex f) {
    if(!f.match(F(w))) {
        cerr << ErrColor << "Not F(w) form: " << f << RESET << endl;
        exit(1);
    }
    if(ccF[f].nops()>0) return ccF[f];
    
    auto al = f.op(0);
    lst sl, rl;
    ex r=0, s=0, sid=0;
    for(int i=0; i<al.nops(); i++) {
        if(al.op(i) > ex(1)/2) {
            sid += pow(2, i);
            r += al.op(i);
            rl.append(al.op(i));
        } else {
            s -= al.op(i);
            sl.append(al.op(i));
        }
    }
    lst ret = lst{sid, r+s, r, s, sl, rl};
    
    if(Cuts.nops()>0) {
        lst clst;
        for(auto c : Cuts) {
            int ci = ex_to<numeric>(c).to_int();
            clst.append(al.op(ci));
        }
        ret.prepend(clst);
    }
    
    ccF[f] = ret;
    ccFinv[ret] = f;
    return ret;
}

lst Laporta::I2lst(ex ibp) {
    if(ccI[ibp].nops()>0) return ccI[ibp];
    
    lst fis;
    exset fiset;
    ibp.find(F(w), fiset);
    for(auto &item : fiset) {
        fis.append(F2lst(item));
    }
    
    for(int i=0; i<fis.nops(); i++) {
        for(int j=i+1; j<fis.nops(); j++) {
            if(ex_less(fis.op(i), fis.op(j))) {
                auto tmp = fis.op(i);
                fis.let_op(i) = fis.op(j);
                fis.let_op(j) = tmp;
            }
        }
    }
    ccImax[ibp] = ex_to<lst>(fis.op(0));
    
    fis.prepend(0);
    fis.let_op(0) = fis.op(1);
    fis.let_op(1) = fis.nops()-1;
    
    ccI[ibp] = fis;
    return fis;
}

ex Laporta::collectF(ex expr) {
    if(!expr.subs(F(w)==0).is_zero()) {
        cerr << ErrColor << "expr is NOT zero with F(w)==0: " << expr << RESET << endl;
        exit(1);
    }
    exset fs;
    expr.find(F(w), fs);
    ex ret = 0;
    for(auto item : fs) {
        ret += normal(expr.coeff(item, 1)) * item;
    }
    return ret;
}

void Laporta::Prepare(lst loop, lst ext, lst prop, lst repl) {

    if(Verbose > 0) cout << "  IBP Prepare @ " << now() << endl;

    lst sps;
    for(auto it : loop) {
        for(auto ii : loop) sps.append(it*ii);
        for(auto ii : ext) sps.append(it*ii);
    }
    sps.sort();
    sps.unique();
    
    lst sp2s, s2sp, ss;
    for(auto item : sps) {
        symbol si;
        ss.append(si);
        sp2s.append(item==si);
        s2sp.append(si==item);
    }
    
    lst eqns;
    for(int i=0; i<prop.nops(); i++) {
        auto eq = prop.op(i).expand();
        eq = eq.subs(sp2s, subs_options::algebraic);
        eq = eq.subs(repl, subs_options::algebraic);
        eqns.append(eq == iWF(i));
    }
    auto s2p = lsolve(eqns, ss);

    preIBPs.remove_all();
    lst ns;
    for(int i=0; i<prop.nops(); i++) ns.append(a(i));
    for(auto l : loop) {
        ex ibp = 0;
        symbol sl;
        for(int i=0; i<prop.nops(); i++) {
            auto ns_tmp = ns;
            ns_tmp.let_op(i) = ns.op(i) + 1;
            auto dp = prop.op(i).subs(l==sl).diff(sl).subs(sl==l);
            ibp -= a(i) * F(ns_tmp) * dp;
        }
        
        for(auto ii : ext) {
            auto ibp_tmp = ibp * ii;
            ibp_tmp = ibp_tmp.expand();
            ibp_tmp = ibp_tmp.subs(sp2s, subs_options::algebraic);
            ibp_tmp = ibp_tmp.subs(repl, subs_options::algebraic);
            ibp_tmp = ibp_tmp.subs(s2p, subs_options::algebraic);
            ex res = 0;
            for(int i=0; i<prop.nops(); i++) {
                auto ci = ibp_tmp.coeff(iWF(i), 1);
                ci = MapFunction([i](const ex & e, MapFunction &self)->ex{
                    if(e.match(F(w))) {
                        auto tmp = e.op(0);
                        tmp.let_op(i) = tmp.op(i)-1;
                        return F(tmp);
                    } else if(!e.has(F(w))) return e;
                    else return e.map(self);
                })(ci);
                res += ci;
            }
            res += ibp_tmp.subs(lst{iWF(w)==0});
            preIBPs.append(res);
        }
        
        for(auto ii : loop) {
            auto ibp_tmp = ibp * ii;
            ibp_tmp = ibp_tmp.expand();
            ibp_tmp = ibp_tmp.subs(sp2s, subs_options::algebraic);
            ibp_tmp = ibp_tmp.subs(repl, subs_options::algebraic);
            ibp_tmp = ibp_tmp.subs(s2p, subs_options::algebraic);
            ex res = 0;
            for(int i=0; i<prop.nops(); i++) {
                auto ci = ibp_tmp.coeff(iWF(i), 1);
                ci = MapFunction([i](const ex &e, MapFunction &self)->ex {
                    if(e.match(F(w))) {
                        auto tmp = e.op(0);
                        tmp.let_op(i) = tmp.op(i)-1;
                        return F(tmp);
                    } else if(!e.has(F(w))) return e;
                    else return e.map(self);
                })(ci);
                res += ci;
            }
            res += ibp_tmp.subs(lst{iWF(w)==0});
            if(ii==l) res += d*F(ns);
            preIBPs.append(res);
        }
    }
    
}

void Laporta::Generate2(lst seed) {
    lst ns, a2s;
    for(int i=0; i<seed.nops(); i++) {
        Symbol s("a"+to_string(i));
        a2s.append(a(i)==s+Symbol("eps"));
        ns.append(s);
    }
    lst seeds;
    for(auto const & item : preIBPs) {
        exset fs;
        item.find(F(w), fs);
        for(auto fi : fs) {
            lst eqs;
            for(int i=0; i<fi.op(0).nops(); i++) {
                eqs.append(subs(fi.op(0).op(i)==seed.op(i),a2s));
            }
            auto sol = lsolve(eqs, ns);
            seeds.append(subs(ns,sol));
            auto ii = item.subs(a2s).subs(sol);
            if(ii.is_zero()) continue;
            IBPs.append(ii);
        }
    }
    
    for(auto seed : seeds)
    for(auto const & item : preIBPs) {
        exset fs;
        item.find(F(w), fs);
        for(auto fi : fs) {
            lst eqs;
            for(int i=0; i<fi.op(0).nops(); i++) {
                eqs.append(subs(fi.op(0).op(i)==seed.op(i),a2s));
            }
            auto sol = lsolve(eqs, ns);
            auto ii = item.subs(a2s).subs(sol);
            if(ii.is_zero()) continue;
            IBPs.append(ii);
        }
    }
    
    if(Cuts.nops()>1) {
        int total = IBPs.nops();
        for(int i=0; i<IBPs.nops(); i++) {
            auto tmp = IBPs.op(i);
            exset fs;
            find(tmp, F(w), fs);
            exmap repl;
            for(auto fi : fs){
                for(auto ic : Cuts) {
                    int j = ex_to<numeric>(ic).to_int();
                    if(fi.op(0).op(j)<=0) {
                        repl[fi]=0;
                        break;
                    }
                }
            }
            IBPs.let_op(i) = tmp.subs(repl);
        }
    }
}

void Laporta::Generate(vector<lst> seeds) {
    if(Verbose > 0) cout << " IBP Generate @ " << now() << endl;

    IBPs.remove_all();
    lst ns;
    for(int i=0; i<seeds[0].nops(); i++) ns.append(a(i));
    for(auto seed : seeds) {
        for(auto const & item : preIBPs) {
            auto ii = item.subs(ns, seed);
            if(ii.is_zero()) continue;
            IBPs.append(ii);
    }}
    
    if(Cuts.nops()>1) {
        int total = IBPs.nops();
        for(int i=0; i<IBPs.nops(); i++) {
            auto tmp = IBPs.op(i);
            exset fs;
            find(tmp, F(w), fs);
            exmap repl;
            for(auto fi : fs){
                for(auto ic : Cuts) {
                    int j = ex_to<numeric>(ic).to_int();
                    if(fi.op(0).op(j)<=0) {
                        repl[fi]=0;
                        break;
                    }
                }
            }
            IBPs.let_op(i) = tmp.subs(repl);
        }
    }
    
    if(Verbose > 0) cout << " - Number of Identities: " << IBPs.nops() << endl << flush;

}

void Laporta::Reduce() {

    if(Verbose > 0) cout << now() << " - IBP Reduce ..." << endl << flush;
    
    if(Verbose > 1) cout << " - Sorting ... @ " << now(false) << endl << flush;
    vector<ex> ibps;
    for(auto item : IBPs) {
        if(item.is_zero()) continue;
        I2lst(item);
        //if(F2lst(F(lst{1,1,1,0,1,1,1})).op(2)+2<ccImax[item].op(2)) continue;
        ibps.push_back(item);
    }

    exset fs;
    find(IBPs, F(w), fs);
    vector<ex> fvec;
    for(auto item : fs) fvec.push_back(item);
        
    // sort Fs
    sort(fvec.begin(), fvec.end(), [this](auto &a, auto &b) {
        return ex_less(F2lst(b),F2lst(a));
    }); // larger first
    
    // sort Identities
    sort(ibps.begin(), ibps.end(), [this](auto &a, auto &b) {
        lst ai = I2lst(a);
        lst bi = I2lst(b);
        return ex_less(bi,ai);
    }); // larger first
    
    map<ex, int, ex_is_less> f2idx;
    for(int i=0; i<fvec.size(); i++) f2idx[fvec[i]] = i;
    
    int cols = fvec.size();
    int rows = ibps.size();
    
    if(Verbose > 1) cout << " - Fermating @ " << now(false) << " :> " << flush;
    
    Fermat fermat;
    fermat.Init();
    
    auto Variables = gather_symbols(ibps);
    ostringstream ss;
    for(auto j : Variables) ss << "&(J="<<j<<");" << endl;
    fermat.Execute(ss.str());
    ss.clear();
    ss.str("");
    ss << "Array mat[" << rows << "," << cols+1 << "] Sparse;" << endl;
    for(int r=0; r<rows; r++) {
        auto ibp = ibps[r];
        exset fs;
        find(ibp, F(w), fs);
        for(auto f : fs) {
            ss << "mat["<<(r+1)<<","<<(f2idx[f]+1)<<"] := " << ibp.coeff(f, 1) << ";" << endl;
            fermat.Execute(ss.str());
            ss.clear();
            ss.str("");
        }
    }

    ss << "Redrowech([mat]);" << endl;
    fermat.Execute(ss.str());
    ss.clear();
    ss.str("");
       
    ss << "&(U=1);" << endl; // ugly printing, the whitespace matters
    ss << "![mat" << endl;
    auto ostr = fermat.Execute(ss.str());
    fermat.Exit();
    
    if(Verbose > 1) cout << now(false) << endl;
    
    // make sure last char is 0
    if(ostr[ostr.length()-1]!='0') throw Error("Reduce: last char is NOT 0.");
    ostr = ostr.substr(0, ostr.length()-1);
    string_trim(ostr);
    
    ostr.erase(0, ostr.find(":=")+2);
    string_replace_all(ostr, "[", "{");
    string_replace_all(ostr, "]", "}");
    Parser fp;
    auto mat2 = fp.Read(ostr);
        
    ex irows[rows];
    for(int i=0; i<rows; i++) {
        irows[i] = 0;
        for(int j=0; j<fvec.size(); j++) {
            if(is_zero(get_op(mat2,i,j))) continue;
            irows[i] += get_op(mat2,i,j) * fvec[j];
        }
    }
    
    if(Verbose > 1) cout << " - Finalizing ... @ " << now(false) << endl << flush;
    
    ibps.clear();
    for(int i=0; i<rows; i++) {
        if(irows[i].is_zero()) continue;
        ibps.push_back(irows[i]);
    }

    Rules.remove_all();
    for(int i=0; i<ibps.size(); i++) {
        auto ibp = ibps[ibps.size()-i-1];
        auto ifmt = I2lst(ibp);
        auto fx = ccFinv[ex_to<lst>(ifmt.op(0))];
        auto c1 = ibp.coeff(fx, 1);
        auto c0 = ibp.coeff(fx, 0);
        auto res = -c0/c1;
        if(!(c1-1).is_zero()) {
            cout << "c1 = " << c1 << endl;
            cout << "ibps.size()-i-1 : " << ibps.size()-i-1 << endl;
            cout << ibp << endl;
            cout << endl << endl;
        }
        Rules.append(fx==res);
    }
    
}

}
