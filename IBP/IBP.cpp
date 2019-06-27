#include "IBP.h"

namespace HepLib {

symbol const IBP::ep("ep");
symbol const IBP::eps("eps");
symbol const IBP::d("d");

lst IBP::formatF(ex f,  FormatCache &cache) {
    assert(f.match(F(wild())));
    if(cache.F2fmt[f].nops()>0) return cache.F2fmt[f];
    
    auto al = f.op(0);
    lst sl, rl;
    ex r=0, s=0, sid=0;
    for(int i=0; i<al.nops(); i++) {
        if(al.op(i).subs(lst{eps==0}) > ex(1)/2) {
            sid += pow(2, i);
            r += al.op(i);
            rl.append(al.op(i));
        } else {
            s -= al.op(i);
            sl.append(al.op(i));
        }
    }
    lst ret = lst{r, s, sl, rl, sid};
    
    if(cache.Cuts.nops()>0) {
        lst clst;
        for(auto c : cache.Cuts) {
            int ci = ex_to<numeric>(c).to_int();
            clst.append(al.op(ci));
        }
        ret.prepend(clst);
    }
    
    cache.F2fmt[f] = ret;
    cache.fmt2F[ret] = f;
    return ret;
}

lst IBP::formatI(ex ibp, FormatCache &cache) {
    if(cache.I2fmt[ibp].nops()>0) return cache.I2fmt[ibp];
    
    lst fis;
    exset fiset;
    ibp.find(F(wild()), fiset);
    for(auto &item : fiset) {
        fis.append(formatF(item, cache));
    }
    
    for(int i=0; i<fis.nops(); i++) {
        for(int j=i+1; j<fis.nops(); j++) {
            if(less(ex_to<lst>(fis.op(i)), ex_to<lst>(fis.op(j)))) {
                auto tmp = fis.op(i);
                fis.let_op(i) = fis.op(j);
                fis.let_op(j) = tmp;
            }
        }
    }
    
    //fis.prepend(0);
    //fis.let_op(0) = fis.op(1);
    //fis.let_op(1) = fis.nops()-1;
    
    cache.I2fmt[ibp] = fis;
    return fis;
}

bool IBP::less(lst nls1, lst nls2) {
    lst ls1 = ex_to<lst>(subs(nls1, lst{eps==0}));
    lst ls2 = ex_to<lst>(subs(nls2, lst{eps==0}));
    if(ls1.nops()<1 || ls2.nops()<1) {
        return ls1.nops() < ls2.nops();
    }
    
    for(int i=0; (i<ls1.nops() && i<ls2.nops()); i++) {
        auto a = ls1.op(i);
        auto b = ls2.op(i);
        if(a.is_equal(b)) continue;
        if(is_a<lst>(a)) {
            assert(is_a<lst>(b));
            return less(ex_to<lst>(a), ex_to<lst>(b));
        } else {
            return a < b;
        }
    }
    
    return ls1.nops() < ls2.nops();
}

ex IBP::collectF(ex expr) {
    assert(expr.subs(F(wild())==0).is_zero());
    exset fs;
    expr.find(F(wild()), fs);
    ex ret = 0;
    for(auto item : fs) {
        ret += normal(expr.coeff(item, 1)) * item;
    }
    return ret;
}

void IBP::Prepare(lst loop, lst ext, lst prop, lst repl) {

    if(Verbose > 0) cout << now() << " - IBP Prepare ..." << endl << flush;

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
        eqns.append(eq == P(i));
    }
    auto s2p = lsolve(eqns, ss);

    preIBPs.remove_all();
    lst ns;
    for(int i=0; i<prop.nops(); i++) ns.append(n(i));
    for(auto l : loop) {
        ex ibp = 0;
        symbol sl;
        for(int i=0; i<prop.nops(); i++) {
            auto ns_tmp = ns;
            ns_tmp.let_op(i) = ns.op(i) + 1;
            auto dp = prop.op(i).subs(l==sl).diff(sl).subs(sl==l);
            ibp -= n(i) * F(ns_tmp) * dp;
        }
        
        for(auto ii : ext) {
            auto ibp_tmp = ibp * ii;
            ibp_tmp = ibp_tmp.expand();
            ibp_tmp = ibp_tmp.subs(sp2s, subs_options::algebraic);
            ibp_tmp = ibp_tmp.subs(repl, subs_options::algebraic);
            ibp_tmp = ibp_tmp.subs(s2p, subs_options::algebraic);
            ex res = 0;
            for(int i=0; i<prop.nops(); i++) {
                auto ci = ibp_tmp.coeff(P(i), 1);
                ci = GiNaC_Replace(ci, F(wild()), [&](auto fi) {
                    auto tmp = fi.op(0);
                    tmp.let_op(i) = tmp.op(i)-1;
                    return F(tmp);
                });
                res += ci;
            }
            res += ibp_tmp.subs(lst{P(wild())==0});
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
                auto ci = ibp_tmp.coeff(P(i), 1);
                ci = GiNaC_Replace(ci, F(wild()), [&](auto fi) {
                    auto tmp = fi.op(0);
                    tmp.let_op(i) = tmp.op(i)-1;
                    return F(tmp);
                });
                res += ci;
            }
            res += ibp_tmp.subs(lst{P(wild())==0});
            if(ii==l) res += d*F(ns);
            preIBPs.append(res);
        }
    }
    
}

void IBP::Generate(vector<lst> seeds) {

    if(Verbose > 0) cout << now() << " - IBP Generate ..." << endl << flush;

    IBPs.remove_all();
    lst ns;
    for(int i=0; i<seeds[0].nops(); i++) ns.append(n(i));
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
            find(tmp, F(wild()), fs);
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

void IBP::Reduce() {

    if(Verbose > 0) cout << now() << " - IBP Reduce ..." << endl << flush;

    Variables.append(d);
    Variables.append(ep);
    Variables.append(eps);
    Variables.sort();
    Variables.unique();
    
    if(Verbose > 1) cout << " - Sorting ... @ " << now(false) << endl << flush;
    vector<ex> ibps;
    for(auto item : IBPs) {
        if(!item.is_zero()) ibps.push_back(item);
    }
    exset fs;
    find(IBPs, F(wild()), fs);
    vector<ex> fvec;
    for(auto item : fs) fvec.push_back(item);
    
    FormatCache cache;
    cache.Cuts = Cuts;
    
    // sort Fs
    sort(fvec.begin(), fvec.end(), [&](auto &a, auto &b) {
        return less(formatF(b, cache), formatF(a, cache));
    }); // larger first
    
    // sort Identities
    sort(ibps.begin(), ibps.end(), [&](auto &a, auto &b) {
        lst aifmt = formatI(a, cache);
        lst bifmt = formatI(b, cache);
        if(less(aifmt, bifmt)) return false;
        if(less(bifmt, aifmt)) return true;
        return !ex_is_less()(a,b);
    }); // larger first
    
    map<ex, int, ex_is_less> f2idx;
    for(int i=0; i<fvec.size(); i++) f2idx[fvec[i]] = i;
    
    int cols = fvec.size();
    int rows = ibps.size();
    
    auto pid = getpid();
    ostringstream foss;
    foss << "fin" << pid;
    const char* ifn = foss.str().c_str();
    foss.clear();
    foss.str("");
    foss << "fout" << pid;
    const char* ofn = foss.str().c_str();
    
    ofstream ofs;
    ofs.open(ifn, ios::out);
    ofs << "&(S="<<ofn<<");" << endl;
    ofs << "&(_d=1000000000);" << endl;
    // &_s: Suppress/don't suppress display of long polynomials.
    ofs << "&_s;" << endl;
    // &_o: http://home.bway.net/lewis/fer64mono.html
    ofs << "&(_o=1000);" << endl;
    // & t: Toggle switch to turn on/off a certain fast probabalistic algorithm
    // to test if one multivariate polynomial divides another over ground ring Z.
    ofs << "&(_t=0);" << endl;
    
    for(auto j : Variables) ofs << "&(J="<<j<<");" << endl;
    ofs << "Array fM[" << rows << "," << cols+1 << "] Sparse;" << endl;
    for(int r=0; r<rows; r++) {
        auto ibp = ibps[r];
        exset fs;
        find(ibp, F(wild()), fs);
        for(auto f : fs) {
            ofs << "fM["<<(r+1)<<","<<(f2idx[f]+1)<<"] := " << ibp.coeff(f, 1) << ";" << endl;
        }
    }
    ofs << "Redrowech([fM]);" << endl;
    // &U: Toggle switch to enable ugly display.
    ofs << "&(U=1);" << endl;
    ofs << "!(&o,[fM]);" << endl;
    ofs << "&q;" << endl;
    ofs << "&x;" << endl;
    ofs.close();
    
    if(Verbose > 1) cout << " - Fermating ... @ " << now(false) << endl << flush;
    RunFermat(ifn);
    
    symtab table;
    for(auto var : Variables) {
        ostringstream oss;
        oss << var;
        table[oss.str()] = var;
    }
    parser reader(table);
    
    ifstream ifs;
    ifs.open(ofn, ios::in);
    
    ex irows[rows];
    for(int i=0; i<rows; i++) irows[i] = 0;
    ostringstream oss;
    int row, col;
    while(true) {
        int c = ifs.get();
        if (c == EOF) break;
        if (isspace(c) || c == '\n') continue;
        if(c == '[') {
            oss.clear();
            oss.str("");
        } else if(c == ',') {
            row = ex_to<numeric>(reader(oss.str())).to_int();
            oss.clear();
            oss.str("");
        } else if(c == ']') {
            col = ex_to<numeric>(reader(oss.str())).to_int();
        } else if(c == '=') {
            oss.clear();
            oss.str("");
        } else if(c == ';') {
            irows[row-1] += reader(oss.str()) * fvec[col-1];
            oss.clear();
            oss.str("");
        } else {
            oss << (char)c;
        }
    }
    ifs.close();
    remove(ifn);
    remove(ofn);
    
    if(Verbose > 1) cout << " - Finalizing ... @ " << now(false) << endl << flush;
    
    ibps.clear();
    for(int i=0; i<rows; i++) {
        if(irows[i].is_zero()) continue;
        ibps.push_back(irows[i]);
    }

    FSolution.remove_all();
    for(int i=0; i<ibps.size(); i++) {
        auto ibp = ibps[ibps.size()-i-1];
        auto ifmt = formatI(ibp, cache);
        auto fx = cache.fmt2F[ex_to<lst>(ifmt.op(0))];
        auto c1 = ibp.coeff(fx, 1);
        auto c0 = ibp.coeff(fx, 0);
        auto res = -c0/c1;
        if(!(c1-1).is_zero()) {
            cout << "c1 = " << c1 << endl;
            cout << "ibps.size()-i-1 : " << ibps.size()-i-1 << endl;
            cout << ibp << endl;
            cout << endl << endl;
        }
        FSolution.append(fx==res);
    }
    
}

REGISTER_FUNCTION(P, dummy())
REGISTER_FUNCTION(F, dummy())
REGISTER_FUNCTION(L, dummy())
REGISTER_FUNCTION(n, dummy())

}
