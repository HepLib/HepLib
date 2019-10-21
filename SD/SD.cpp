#include "SD.h"
#include <math.h>

namespace HepLib {

symbol const SD::ep("ep");
symbol const SD::eps("eps");
symbol const SD::vs("s");
symbol const SD::vz("z");
symbol const SD::iEpsilon("iEpsilon");
realsymbol const SD::NaN("NaN");
bool SD::use_dlclose = true;
bool SD::debug = false;

vector<exmap> SecDecBase::x2y(const lst &xpols) {
    ex xpol = 1;
    for(int i=0; i<xpols.nops(); i++) xpol *= xpols.op(i);
    return x2y(xpol);
}

/*-----------------------------------------------------*/
/*                    CheckFAtX1                         */
/*-----------------------------------------------------*/
bool SD::IsBadF1(ex f, vector<exmap> vmap) {
    for(auto &vi : vmap) {
        auto ft = f.subs(vi);
        auto xs_tmp = get_x_from(ft);
        auto ys_tmp = get_y_from(ft);
        int ysn = ys_tmp.size();
        for(int i=0; i<xs_tmp.size(); i++) {
            vi[xs_tmp[i]] = y(ysn+i);
        }

        // need collect_common_factors
        ft = collect_common_factors(ft.expand());

        if(is_exactly_a<mul>(ft)) {
            ex ret = 1;
            for (auto item : ft) {
                if( !(item.match(y(wild())) || item.match(pow(y(wild()), wild(1)))) ) {
                    ret *= item;
                }
            }
            ft = ret;
        } else if ( ft.match(y(wild())) || ft.match(pow(y(wild()), wild(1))) ) {
            ft = 1;
        }

        ys_tmp = get_y_from(ft);
        for(int pi=1; pi<std::pow(2, ys_tmp.size()); pi++) {
            lst yrepl;
            int ci = pi;
            for(int i=0; i<ys_tmp.size(); i++) {
                yrepl.append(ys_tmp[i] == ((ci%2)==1 ? 1 : 0));
                ci /= 2;
            }
            if(normal(ft.subs(yrepl)).is_zero()) {
                return true;
            }
        }
    }
    return false;
}

vector<pair<lst, lst>> SD::AutoF1(pair<lst, lst> po_ex) {
    assert(Deltas.size()<1);
    lst const exlist = po_ex.second;
    assert((exlist.op(0)-1).is_zero());
    auto xs = get_x_from(po_ex.first);
    if(xs.size()<1) xs = get_y_from(po_ex.first);
    int nx = xs.size();
    for(int nn=0; nn<=nx; nn++) {
    for(int pi=0; pi<std::pow(2, nx); pi++) {
        int cpi = pi, cn1 = 0;
        for(int i=0; i<nx; i++) {
            if((cpi % 2) == 1) cn1++;
            cpi /= 2;
        }
        if(cn1 != nn) continue;
        
        lst polists = lst{ po_ex.first };
        int bi = 0, bs = BisectionPoints.nops();
        cpi = pi;
        for(int i=0; i<nx; i++) {
            if((cpi % 2) == 1) {
                lst polists2;
                auto x0 = BisectionPoints.op(bi % bs);
                for(auto item : polists) {
                    symbol xy;
                    auto tmp = ex_to<lst>(subs(subs(item, lst{xs[i]==x0*xy}), lst{xy==xs[i]}));
                    tmp.let_op(0) = tmp.op(0) * x0;
                    polists2.append(tmp);
                    
                    tmp = ex_to<lst>(subs(subs(item, lst{xs[i]==(x0-1)*xy+1}), lst{xy==xs[i]}));
                    tmp.let_op(0) = tmp.op(0) * (1-x0);
                    polists2.append(tmp);
                }
                polists = polists2;
                bi++;
            }
            cpi /= 2;
        }
        bool OK = true;
        Digits = 50;
        for(auto polist : polists) {
            lst sdList;
            for(int i=0; i<polist.nops(); i++) {
                auto tmp = polist.op(i);
                auto ntmp = exlist.op(i);
                if(!tmp.subs(lst{x(wild())==0, y(wild())==0}).normal().is_zero()) continue;
                if( (!tmp.has(x(wild())) && !tmp.has(y(wild()))) || (is_a<numeric>(ntmp) && ntmp.evalf()>0) ) continue;
                sdList.append(tmp);
            }
            vector<exmap> vmap = SecDec->x2y(sdList);
            if(IsBadF1(polist.op(1), vmap)) {
                OK = false;
                break;
            }
        }
        
        if(OK) {
            vector<pair<lst, lst>> res;
            for(auto item : polists) res.push_back(make_pair(ex_to<lst>(item), exlist));
            return res;
        }
    }}
    
    cerr << RED << "F: " << po_ex.first.op(1) << RESET << endl;
    cerr << RED << "F1 Failed with ALL possible bisections!" << RESET << endl;
    assert(false);
    return vector<pair<lst, lst>>();
}

/*-----------------------------------------------------*/
/*                       SD                            */
/*-----------------------------------------------------*/
//return a lst, element pattern: { {{x1,n1}, {x2,n2}, ...}, {{e1, n1},{e2,n2}, ...} }
//e1 is a const term, e2 still the F-term
vector<lst> SD::DS(pair<lst, lst> po_ex) {
    // 1st element in input polist is the constant term, guess NOT necessary
    // 2nd element in input polist is the F-term, required!
    lst const polist = po_ex.first;
    lst const exlist = po_ex.second;
    lst sdList;
    for(int i=0; i<polist.nops(); i++) {
        auto tmp = polist.op(i);
        auto ntmp = exlist.op(i);
        if(!tmp.subs(lst{x(wild())==0, y(wild())==0}).normal().is_zero()) continue;
        Digits = 50;
        if( (!tmp.has(x(wild())) && !tmp.has(y(wild()))) || (is_a<numeric>(ntmp) && ntmp.evalf()>0) ) continue;
        sdList.append(tmp);
    }

    vector<exmap> vmap = SecDec->x2y(sdList);

    vector<lst> sd;
    for(auto &vi : vmap) {
        auto ypolist = polist.subs(vi);
        auto xs_tmp = get_x_from(ypolist);
        auto ys_tmp = get_y_from(ypolist);
        int ysn = ys_tmp.size();
        for(int i=0; i<xs_tmp.size(); i++) {
            vi[xs_tmp[i]] = y(ysn+i);
        }
        if(xs_tmp.size()>0) ypolist = polist.subs(vi);

        // need collect_common_factors
        auto ft = collect_common_factors(ypolist.op(1).expand());
        
        ex ct = 1, fsgin = 1;
        if(is_a<mul>(ft)) {
            ex ret = 1;
            for (auto item : ft) {
                if( !(item.match(y(wild())) || item.match(pow(y(wild()), wild(1)))) ) {
                    ret *= item;
                }
            }
            ft = ret;
        } else if ( ft.match(y(wild())) || ft.match(pow(y(wild()), wild(1))) ) {
            ft = 1;
        }
        
        bool need_contour_deformation = ft.has(PL(wild()));
        if(ft.has(y(wild())) && !need_contour_deformation) {
            auto tmp = ft.subs(nReplacements).expand();
            if(is_a<add>(tmp)) {
                need_contour_deformation = false;
                auto first = tmp.op(0).subs(y(wild())==1);
                for(auto item : tmp) {
                    assert(is_a<numeric>(item.subs(y(wild())==1)*first));
                    if(item.subs(y(wild())==1)*first<0) {
                        need_contour_deformation = true;
                        break;
                    }
                }
            }
        }

        if(!need_contour_deformation) {
            auto tmp = ft.subs(y(wild())==1).subs(nReplacements);
            if(!is_a<numeric>(tmp)) {
                cerr << "tmp = " << tmp << endl;
                cerr << "tmp is NOT a numeric." << endl;
                assert(false);
            }
            Digits = 50;
            if( tmp.evalf() < 0 ) {
                ct = exp(-I * Pi * exlist.op(1));
                fsgin = -1;
            }
            ft = 1;
        } else {
            auto tmp = ft.subs(nReplacements);
            auto tn = tmp.subs(y(wild())==numeric("1/17"));
            for(auto kv : ParameterUB) {
                int k = kv.first;
                tn = tn.subs(PL(k)==(ParameterUB[k]+ParameterLB[k])/ex(2));
            }
            if(tn.has(PL(wild()))) {
                cerr << "PL still exists in " << tn << endl;
                assert(false);
            }
            Digits = 50;
            if(!is_a<numeric>(tn.evalf())) {
                cerr << "tn is not numeric: " << tn.eval() << endl;
                assert(false);
            }
            
            if(tn.evalf() < 0) tmp = ex(0)-tmp;
            double tmin = FindMinimum(tmp, true);
            if(tmin > 1E-5) {
                if(tn.evalf() < 0) {
                    ct = exp(-I * Pi * exlist.op(1));
                    fsgin = -1;
                }
                ft = 1;
                need_contour_deformation = false;
            } 
        }

        exmap ymol;
        
        // need collect_common_factors
        auto det = collect_common_factors(vi[x(-1)]);
        assert(!is_exactly_a<add>(det));
        auto ys = get_xy_from(det);
        ex det1 = 1;
        for(int i=0; i<ys.size(); i++) {
            ymol[ys[i]] = ymol[ys[i]] + det.degree(ys[i]);
            det1 *= pow(ys[i], det.degree(ys[i]));
        }
        assert(is_a<numeric>(det/det1) && ex_to<numeric>(det/det1).is_integer());
        
        lst ftxlst = lst{0};
        for(auto xi : get_xy_from(ft)) ftxlst.append(xi);
        
        lst pol_exp_lst;
        pol_exp_lst.append(lst{ subs(FTX(ft,ftxlst)*CT(ct*det/det1), y(wild())==x(wild())), 1 });

        for(int i=0; i<ypolist.nops(); i++) {
            auto tmp = ypolist.op(i);
            auto nex = exlist.op(i);
            bool nchk = (!is_a<numeric>(nex) || ex_to<numeric>(nex)<0);
    
            if(i==1 && fsgin<0) tmp = -tmp; // check complete negtive F-term
            auto tmp1 = tmp;
            while(true) {
                tmp=tmp1.subs(pow(y(wild(1))*wild(2), wild(3))==pow(y(wild(1)),wild(3)) * pow(wild(2), wild(3)));
                if((tmp1-tmp).is_zero()) break;
                tmp1 = tmp;
            }
            
            // need collect_common_factors
            if(tmp.has(y(wild()))) tmp = collect_common_factors(tmp.expand());

            lst tmps;
            if(is_exactly_a<mul>(tmp)) {
                for (auto item : tmp) tmps.append(item);
            } else {
                tmps.append(tmp);
            }
            
            ex rem = 1;
            ex ct = 1;
            for (auto item : tmps) {
                if( item.match(y(wild())) || item.match(pow(y(wild()), wild(1))) ) {
                    auto yi = get_xy_from(item)[0];
                    ymol[yi] = ymol[yi] + (item.nops()<2 ? 1 : item.op(1)) * exlist.op(i);
                } else if(!item.has(y(wild())) && !item.has(x(wild()))) {
                    if(is_a<numeric>(exlist.op(i)) && ex_to<numeric>(exlist.op(i)).is_integer()) {
                        ct *= item;
                    } else if(!item.has(PL(wild()))) {
                        auto tr = item.subs(nReplacements);
                        Digits = 50;
                        if(!is_a<numeric>(tr.evalf())) {
                            cerr << "not numeric - item: " << tr << " ; " << item << endl;
                            assert(false);
                        }
                        auto nitem = ex_to<numeric>(tr.evalf());
                        if( nitem.is_real() && nitem<0 ) {
                            ct *= -ex(1)*item;
                            rem *= -ex(1);
                        } else {
                            ct *= item;
                        }
                    } else {
                        rem *= item;
                    }
                } else {
                    if(nchk && item.subs(y(wild())==0).subs(iEpsilon==0).normal().is_zero()) {
                        cerr << "zero - item: " << item << endl;
                        cerr << "exlist.op(i) = " << exlist.op(i) << endl;
                        assert(false);
                    }
                    rem *= item;
                }
            }
            
            ex pnp = rem-(i==1 && ft!=1 ? iEpsilon : ex(0));
            pnp = pnp.subs(y(wild())==x(wild()));
            ex pnn = exlist.op(i);
            pnn = pnn.subs(y(wild())==x(wild()));
            
            pol_exp_lst.append(lst{pnp, pnn});
            pol_exp_lst.let_op(0) = pol_exp_lst.op(0).subs(CT(wild()) == CT(wild()*pow(ct, exlist.op(i)))).subs(CT(0)==0);
        }

        lst x_n_lst;
        for(auto & kv : ymol) {
            auto k = kv.first.subs(y(wild())==x(wild()));
            auto v = kv.second.subs(y(wild())==x(wild()));
            if(is_a<numeric>(v)) {
                auto nv = ex_to<numeric>(v);
                if(nv<=-1) {
                    cerr << endl << "nv: " << nv << ", NOT larger than -1." << endl;
                    assert(false);
                }
            }
            x_n_lst.append(lst{k, v});
        }
        x_n_lst.sort();
        sd.push_back(lst{x_n_lst, pol_exp_lst});
    }
    
    return sd;
}

// 1st element in returned lst1 is the constant term
// 2nd element in both returned and inputed lst1 is the F-term
pair<lst, lst> SD::Normalize(const pair<lst, lst> &input) {
    ex const_term = 1;
    lst plst, nlst;
    for(int i=0; i<input.first.nops(); i++) {
        if(input.second[i].is_zero() || input.first[i]==ex(1)) continue;
        if(i!=1 && !input.first[i].has(x(wild())) && !input.first[i].has(y(wild()))) {
            const_term *= pow(input.first[i], input.second[i]);
        } else {
            auto ptmp = input.first[i];
            auto ntmp = input.second[i];
            if(is_exactly_a<mul>(ptmp)) {
                ex tmul = 1;
                for(int j=0; j<ptmp.nops(); j++) {
                    auto tmp = ptmp.op(j);
                    if(!tmp.has(x(wild())) && !tmp.has(y(wild()))) {
                        if(is_a<numeric>(ntmp) && ex_to<numeric>(ntmp).is_integer()) {
                            const_term *=  pow(tmp,ntmp);
                        } else if((tmp-vs).is_zero() || tmp.match(pow(vs,wild()))) {
                            const_term *=  pow(tmp,ntmp);
                        } else if(!tmp.has(PL(wild()))) {
                            auto tr = tmp.subs(nReplacements);
                            if(!is_a<numeric>(tr)) {
                                cerr << "tmp: " << tmp << endl;
                                cerr << "nReplacements: " << nReplacements << endl;
                                cerr << "tmp is NOT numeric with nReplacements." << endl;
                                assert(false);
                            }
                            if(ex_to<numeric>(tr)>0) {
                                const_term *=  pow(tmp,ntmp);
                            } else {
                                const_term *= pow(-tmp,ntmp);
                                tmul *= -1;
                            }
                        } else {
                            tmul *= tmp;
                        }
                    } else {
                        tmul *= tmp;
                    }
                }
                if(tmul != 1) {
                    if(i==1) {
                        plst.prepend(tmul);
                        nlst.prepend(ntmp);
                    } else {
                        plst.append(tmul);
                        nlst.append(ntmp);
                    }
                }
            } else {
                if(i == 1) {
                    plst.prepend(ptmp);
                    nlst.prepend(ntmp);
                } else {
                    plst.append(ptmp);
                    nlst.append(ntmp);
                }
            }
        }
    }
    
    plst.prepend(const_term);
    nlst.prepend(1);
    
    lst plst_comb, nlst_comb;
    exmap np;
    map<ex,int,ex_is_less> inp;
    for(int i=0; i<nlst.nops(); i++) np[plst[i]] += nlst[i];
    for(int i=0; i<nlst.nops(); i++) {
        if(inp[plst[i]] != 0) continue;
        plst_comb.append(plst[i]);
        nlst_comb.append(np[plst[i]]);
        inp[plst[i]] = 1;
    }
    
    return make_pair(plst_comb, nlst_comb);
}

void Replacements2(exmap &repl) {
    auto tmp = repl;
    for(auto &kv : repl) {
        kv.second = kv.second.subs(tmp, subs_options::algebraic);
    }
}

bool xPositive(ex expr) {
    auto tmp = expr.expand();
    if(tmp.is_zero()) return true;
    bool ret = false;
    if(is_a<add>(tmp)) {
        for(auto item : tmp) {
            auto nit = item.subs(x(wild())==1).normal();
            if(!(is_a<numeric>(nit) && ex_to<numeric>(nit).is_positive())) {
                return false;
            }
        }
        ret = true;
    } else {
        auto ntmp = tmp.subs(x(wild())==1).normal();
        ret = (is_a<numeric>(ntmp) && ex_to<numeric>(ntmp).is_positive());
    }
    return ret;
}

void SD::Initialize(FeynmanParameter fp) {
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    if(fp.Propagators.nops() != fp.Exponents.nops()) {
        cerr << "the length of Propagators and Exponents are NOT equal." << endl;
        assert(false);
    }
    
    if(fp.Prefactor.is_zero()) {
        IsZero = true;
        return;
    }
    
    IsZero = false;
    
    for(auto kv: fp.lReplacements) assert(!(lst{kv.first, kv.second}).has(iEpsilon));
    for(auto kv: fp.tReplacements) assert(!(lst{kv.first, kv.second}).has(iEpsilon));
    for(auto kv: fp.nReplacements) assert(!(lst{kv.first, kv.second}).has(iEpsilon));
    
    auto sop = subs_options::algebraic;
    
    auto ps = fp.Propagators;
    auto ns = fp.Exponents;
    
    auto ls = fp.LoopMomenta;
    auto tls = fp.tLoopMomenta;
    
    Replacements2(fp.lReplacements);
    Replacements2(fp.tReplacements);
    Replacements2(fp.nReplacements);
    
    auto lsubs = fp.lReplacements;
    auto tsubs = fp.tReplacements;
    auto nsubs = fp.nReplacements;
    lst plsubs;
    for(auto kv : ParameterUB) {
        int k = kv.first;
        plsubs.append(PL(k)==(ParameterUB[k]+ParameterLB[k])/ex(2));
    }
    nReplacements = fp.nReplacements;
    
    if(Verbose > 0) cout << now() << " - Initialize ..." << endl << flush;
    
    ex asgn = 1;
    ex a = 0;
    ex ad = (4-2*ep)*ls.nops() + (2-2*ep)*tls.nops();
    int xn = ps.nops();
    ex rem = 0;
    exmap xtNeg;
    
    ex pre = fp.Prefactor; // come from below
    for(int i=0; i<ps.nops(); i++) {
        bool ltQ = false; {
            auto tps = ps.op(i).expand().subs(lsubs,sop).subs(tsubs,sop);
            for(auto lsi : ls) {
                if(tps.has(lsi)) {
                    ltQ = true;
                    break;
                }
            }
            
            if(!ltQ) {
                for(auto lsi : tls) {
                    if(tps.has(lsi)) {
                        ltQ = true;
                        break;
                    }
                }
            }
        }
        
        if(ltQ) a += ns.op(i);
        ex sgn = 0;
        
        if(!ltQ) {
            pre = pre * pow(ps.op(i).expand().subs(lsubs,sop).subs(tsubs,sop), ex(0)-ns.op(i));
            ns.let_op(i) = 0;
            ps.let_op(i) = 1;
            continue;
        } else if(ns.op(i) <= 0) {
            xtNeg[x(i)]=0;
            if(ns.op(i)<0) {
                asgn *= pow(-1, ns.op(i));
                rem += x(i) * ps.op(i).subs(iEpsilon==0);
            }
            continue;
        }

        auto p = ps.op(i).expand().subs(lsubs,sop).subs(tsubs,sop).subs(nsubs).subs(plsubs);
        p = p.subs(lsubs,sop).subs(tsubs,sop).subs(nsubs).subs(plsubs);
        // check loop^2
        for(auto m : ls) {
            assert(is_a<numeric>(p.coeff(m,2)));
            numeric nm = ex_to<numeric>(p.coeff(m,2));
            if(nm.is_zero()) continue;
            sgn = nm>0 ? -1 : 1;
            break;
        }
        // check iEpsilon
        if(sgn.is_zero()) {
            assert(is_a<numeric>(p.coeff(iEpsilon)));
            numeric nm = ex_to<numeric>(p.coeff(iEpsilon));
            if(!nm.is_zero()) sgn = nm>0 ? -1 : 1;
        }
        // check tloop^2
        if(sgn.is_zero()) {
            for(auto m : tls) {
                if(!is_a<numeric>(p.coeff(m,2))) {
                    cerr << "not numeric: " << p.coeff(m,2) << endl;
                    assert(false);
                }
                numeric nm = ex_to<numeric>(p.coeff(m,2));
                if(nm.is_zero()) continue;
                sgn = nm>0 ? -1 : 1;
                break;
            }
        }
        // others
        if(sgn.is_zero()) {
            sgn = 1;
            if(is_a<numeric>(p) && ex_to<numeric>(p)>0) sgn = -1;
            cout << " - Warning: Can NOT determine the iEpsilon sign." << endl;
            cout << " - " << p << " from " << ps.op(i) << endl;
        }
        
        p = (ps.op(i)*sgn).subs(iEpsilon==0);
        if(sgn==-1) asgn *= exp(I * Pi * ns.op(i));
        rem += x(i) * p;
    }

    rem = rem.expand();
    lst uList1, uList2;
    ex u=1, cu=1;
    
    // Loop
    if(ls.nops()>0) {
        u=1;
        for(int i=0; i<ls.nops(); i++) {
            auto t2 = rem.coeff(ls.op(i),2);
            auto t1 = rem.coeff(ls.op(i),1);
            auto t0 = rem.coeff(ls.op(i),0);
            u *= (-t2);
            if(t2==0) {
                IsZero = true;
                return;
            }
            rem = expand(t0 - pow(t1,2)/(4*t2));
        }
        rem = normal(rem.subs(lsubs,sop).subs(lsubs,sop));
        u = normal(u.subs(lsubs,sop));
        for(auto m: tls) assert(!u.has(m));
        
        cu *= u;
        auto u_nd = numer_denom(u);
        ex usgn = u_nd.op(1).subs(xtNeg).subs(x(wild())==ex(1)/2).subs(nsubs).subs(plsubs);
        if(usgn.is_zero()) usgn = u_nd.op(1).subs(xtNeg).subs(x(wild())==ex(1)/3).subs(nsubs).subs(plsubs);
        assert(!usgn.is_zero());
        usgn = normal(usgn)>0 ? 1 : -1;
        
        if(!xPositive(normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacements).subs(plsubs))) {
            cerr << "NOT positive - un: " << normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacements).subs(plsubs) << endl;
            assert(false);
        }
        if(!xPositive(normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacements).subs(plsubs))) {
            cerr << "NOT positive - ud: " << normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacements).subs(plsubs) << endl;
            assert(false);
        }
        
        uList1.append(usgn*u_nd.op(0));
        uList2.append(-(4-2*ep)/2);
        if((usgn*u_nd.op(0)) != 1) {
            uList1.append(usgn*u_nd.op(1));
            uList2.append((4-2*ep)/2);
        }
    } else {
        rem = normal(rem.subs(lsubs,sop).subs(lsubs,sop));
    }
    
    // t-Loop
    if(tls.nops()>0) {
        u=1;
        for(int i=0; i<tls.nops(); i++) {
            auto t2 = rem.coeff(tls.op(i),2);
            auto t1 = rem.coeff(tls.op(i),1);
            auto t0 = rem.coeff(tls.op(i),0);
            u *= (-t2);
            if(t2.is_zero()) {
                IsZero = true;
                return;
            }
            rem = expand(t0 - pow(t1,2)/(4*t2));
        }
        rem = normal(rem.subs(tsubs,sop));
        u = normal(u.subs(lsubs,sop));
        for(auto m: tls) assert(!u.has(m));
        
        cu *= u;
        auto u_nd = numer_denom(u);
        ex usgn = u_nd.op(1).subs(xtNeg).subs(x(wild())==ex(1)/2).subs(nsubs).subs(plsubs);
        if(usgn.is_zero()) usgn = u_nd.op(1).subs(xtNeg).subs(x(wild())==ex(1)/3).subs(nsubs).subs(plsubs);
        assert(!usgn.is_zero());
        usgn = normal(usgn)>0 ? 1 : -1;
        
        if(!xPositive(normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacements).subs(plsubs))) {
            cerr << "NOT positive - un: " << normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacements).subs(plsubs) << endl;
            assert(false);
        }
        if(!xPositive(normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacements).subs(plsubs))) {
            cerr << "NOT positive - ud: " << normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacements).subs(plsubs) << endl;
            assert(false);
        }
        
        uList1.append(usgn*u_nd.op(0));
        uList2.append(-(2-2*ep)/2);
        if(usgn*u_nd.op(1) != 1) {
            uList1.append(usgn*u_nd.op(1));
            uList2.append((2-2*ep)/2);
        }
    } 

    u = normal(cu);
    auto u_nd = numer_denom(u);
    rem = normal(rem * u);
    auto rem_nd = numer_denom(rem);
    
    ex usgn = u_nd.op(1).subs(xtNeg).subs(x(wild())==ex(1)/2).subs(nsubs).subs(plsubs);
    usgn = normal(usgn)>0 ? 1 : -1;
    ex fsgn = rem_nd.op(1).subs(xtNeg).subs(x(wild())==ex(1)/2).subs(nsubs).subs(plsubs);
    fsgn = normal(fsgn)>0 ? 1 : -1;
    
    lst fList1, fList2;
    fList1.append(usgn*u_nd.op(0));
    fList2.append(a-ad/2);
    fList1.append(fsgn*rem_nd.op(0));
    fList2.append(-a+ad/2);
    if(usgn*u_nd.op(1) != 1) {
        fList1.append(usgn*u_nd.op(1));
        fList2.append(-a+ad/2);
    }
    if(fsgn*rem_nd.op(1) != 1) {
        fList1.append(fsgn*rem_nd.op(1));
        fList2.append(a-ad/2);
    }
    
    for(int i=0; i<uList1.nops(); i++) {
        fList1.append(uList1[i]);
        fList2.append(uList2[i]);
    }
    
    vector<pair<lst, lst>> ret;
    ret.push_back(make_pair(fList1, fList2));

    // negative index
    for(int i=0; i<xn; i++) {
    if(is_a<numeric>(ns[i]) && ns[i]<0) {
        for(int j=0; j<-ns[i]; j++) {
            vector<pair<lst, lst>> nret;
            for(auto kv : ret) {
                auto plst = kv.first;
                auto nlst = kv.second;
                for(int ij=0; ij<nlst.nops(); ij++) {
                    auto dtmp = nlst.op(ij) * mma_diff(plst.op(ij),x(i),1,false);
                    if(dtmp.is_zero()) continue;
                    auto plst2 = plst;
                    auto nlst2 = nlst;
                    if((nlst.op(ij)-1).is_zero()) {
                        plst2.let_op(ij) = dtmp;
                    } else {
                        nlst2.let_op(ij) = nlst.op(ij)-1;
                        int nn = plst.nops();
                        if(!(nlst.op(nn-1)-1).is_zero()) {
                            plst2.append(dtmp);
                            nlst2.append(1);
                        } else plst2.let_op(nn-1) = plst.op(nn-1) * dtmp;
                    }
                    nret.push_back(make_pair(plst2, nlst2));
                }
            }
            ret = nret;
        }
        
        for(auto &kv : ret) {
            lstHelper::map_inplace(kv.first, [&](auto &&e) { return e.subs(x(i)==0); });
            lstHelper::map_inplace(kv.second, [&](auto &&e) { return e.subs(x(i)==0); });
        }
    }}

    // simplification
    // ex pre = fp.Prefactor; // moved to above
    pre *= asgn * pow(I,ls.nops()) * pow(Pi, ad/2) * tgamma(a-ad/2);
    for(int i=0; i<ns.nops(); i++) {
        if(is_a<numeric>(ns[i]) && ns.op(i)<=0) continue;
        pre /= tgamma(ns.op(i));
    }
    if(tls.nops()>0) pre *= exp(I * Pi * tls.nops()*(2-2*ep)/2);
    
    ex xpre = 1;
    for(int i=0; i<ns.nops(); i++) {
        if(is_a<numeric>(ns[i]) && ns.op(i)<=1) continue;
        if(is_a<numeric>(ns.op(i))) xpre *= pow(x(i), ns.op(i)-1);
        else for(auto &kv : ret) {
                kv.first.append(x(i));
                kv.second.append(ns.op(i)-1);
             }
    }

    for(auto &kv : ret) {
        kv.first.append(pre);
        kv.second.append(1);
        if(xpre != 1) {
            kv.first.append(xpre);
            kv.second.append(1);
        }
        lstHelper::map_inplace(kv.first, [&](auto &&e) { return collect_common_factors(e); });
        lstHelper::map_inplace(kv.second, [&](auto &&e) { return collect_common_factors(e); });
    }

    lst delta;
    for(int i=0; i<ns.nops(); i++) {
        if(is_a<numeric>(ns[i]) && ns.op(i)<=0) continue;
        delta.append(x(i));
    }
    Deltas.push_back(delta);
    FunExp = ret;

    // Do Other Simplifications
    Normalizes();
    XReOrders();
    Normalizes();
}

/*-----------------------------------------------------*/
/*               's Funtions in SD                       */
/*-----------------------------------------------------*/
void SD::XReOrders() {
if(IsZero) return;
if(Deltas.size()>0) {
    
    exmap xmap;
    for(auto kv : FunExp) {
        lst xExplst;
        for(auto kvf : kv.first) xExplst.append(kvf);
        for(auto kvs : kv.second) xExplst.append(kvs);
        exset xset;
        for(int i=0; i<xExplst.nops(); i++) {
            auto pol = xExplst.op(i);
            pol.find(x(wild()), xset);
            for(auto it=xset.begin(); it!=xset.end(); it++) xmap[*it]++;
        }
    }
    for(auto delta : Deltas) {
        lst xExplst;
        for(auto di : delta) xExplst.append(di);
        exset xset;
        for(int i=0; i<xExplst.nops(); i++) {
            auto pol = xExplst.op(i);
            pol.find(x(wild()), xset);
            for(auto it=xset.begin(); it!=xset.end(); it++) xmap[*it]++;
        }
    }

    vector<ex> xs;
    for(auto kv : xmap) xs.push_back(kv.first);
    sort(xs.begin(), xs.end(), [&](const auto &a, const auto &b){
        return ex_to<numeric>(normal((b-a)).subs(lst{x(wild())==wild()})).is_positive();
    });
    
    lst x2y;
    for(int i=0; i<xs.size(); i++) {
        x2y.append(xs[i]==y(i));
    }
    
    for(auto &kv : FunExp) {
        lstHelper::map_inplace(kv.first, [&](auto &&e) { return e.subs(x2y).subs(y(wild())==x(wild())); });
        lstHelper::map_inplace(kv.second, [&](auto &&e) { return e.subs(x2y).subs(y(wild())==x(wild())); });
    }
    
    for(auto &delta : Deltas) {
        lstHelper::map_inplace(delta, [&](auto &&e) { return e.subs(x2y).subs(y(wild())==x(wild())); });
    }
} else if(Integrands.size()<1) {
    for(auto &kv : FunExp) {
        lst xExplst;
        for(auto kvf : kv.first) xExplst.append(kvf);
        for(auto kvs : kv.second) xExplst.append(kvs);
        exmap xmap;
        exset xset;
        for(int i=0; i<xExplst.nops(); i++) {
            auto pol = xExplst.op(i);
            pol.find(x(wild()), xset);
            for(auto it=xset.begin(); it!=xset.end(); it++) xmap[*it]++;
        }
        
        vector<ex> xs;
        for(auto kv : xmap) xs.push_back(kv.first);
        sort(xs.begin(), xs.end(), [&](const auto &a, const auto &b){
            return ex_to<numeric>(normal((b-a)).subs(lst{x(wild())==wild()})).is_positive();
        });
        
        lst x2y;
        for(int i=0; i<xs.size(); i++) {
            x2y.append(xs[i]==y(i));
        }
        
        lstHelper::map_inplace(kv.first, [&](auto &&e) { return e.subs(x2y).subs(y(wild())==x(wild())); });
        lstHelper::map_inplace(kv.second, [&](auto &&e) { return e.subs(x2y).subs(y(wild())==x(wild())); });
    
    }
} else {
    for(auto & vint : Integrands) {
        exset xset;
        vint.find(x(wild()), xset);
        map<ex, int, ex_is_less> xmap;
        for(auto it=xset.begin(); it!=xset.end(); it++) xmap[*it]++;
        
        vector<ex> xs;
        for(auto kv : xmap) xs.push_back(kv.first);
        sort(xs.begin(), xs.end(), [&](const auto &a, const auto &b){
            return ex_to<numeric>(normal((b-a)).subs(lst{x(wild())==wild()})).is_positive();
        });
        
        lst x2y;
        for(int i=0; i<xs.size(); i++) {
            x2y.append(xs[i]==y(i));
        }
        
        vint = vint.subs(x2y).subs(y(wild())==x(wild()));
    }
}
}

void SD::Normalizes() {
    if(IsZero) return;

    vector<pair<lst, lst>> funexp;
    for(auto fe : FunExp) {
        funexp.push_back(Normalize(fe));
    }
    FunExp.clear();
    FunExp.shrink_to_fit();
    
    exmap fn, ifn;
    for(auto fe : funexp) {
        ex key = 1;
        for(int i=1; i<fe.first.nops(); i++) key *= pow(fe.first, fe.second);
        fn[key] += fe.first.op(0);
    }
    
    for(auto fe : funexp) {
        ex key = 1;
        for(int i=1; i<fe.first.nops(); i++) key *= pow(fe.first, fe.second);
        if(ifn[key]>0) continue;
        auto kv = make_pair(lstHelper::sub(fe.first,1,-1),lstHelper::sub(fe.second,1,-1));
        kv.first.prepend(fn[key]);
        kv.second.prepend(1);
        FunExp.push_back(kv);
        ifn[key] = 1;
    }
}

void SD::XTogethers() {

    vector<pair<lst, lst>> funexp;
    for(auto fe : FunExp) {
        funexp.push_back(fe);
    }
    FunExp.clear();
    FunExp.shrink_to_fit();
    
    map<ex, ex, ex_is_less> fe_cc;
    for(auto fe : funexp) {
        lst key;
        key.append(lst {fe.first.op(1), fe.second.op(1)});
        ex rem = 1;
        for(int i=0; i<fe.second.nops(); i++) {
            if(i==1) continue;
            if(!fe.second.op(i).has(ep) && !fe.second.op(i).has(eps)) {
                if(is_a<numeric>(fe.second.op(i)) && ex_to<numeric>(fe.second.op(i)).is_nonneg_integer()) {
                    rem *= pow(fe.first.op(i), fe.second.op(i));
                    continue;
                }
            }
            key.append(lst{ fe.first.op(i), fe.second.op(i) });
        }
        fe_cc[key] += rem;
    }
    
    for(auto kv : fe_cc) {
        lst klst = lst{1};
        lst vlst = lst{1};
        for(auto fe : kv.first) {
            klst.append(fe.op(0));
            vlst.append(fe.op(1));
        }
        klst.append(kv.second);
        vlst.append(1);

        FunExp.push_back(make_pair(klst, vlst));
    }
    
}

void SD::XExpands() {

    vector<pair<lst, lst>> funexp;
    for(auto fe : FunExp) {
        funexp.push_back(fe);
    }
    FunExp.clear();
    FunExp.shrink_to_fit();
    
    for(auto fe : funexp) {
        lst klst = lst{fe.first.op(0), fe.first.op(1)};
        lst vlst = lst{fe.second.op(0), fe.second.op(1)};
        ex rem = 1;
        for(int i=2; i<fe.second.nops(); i++) {
            if(!fe.second.op(i).has(ep) && !fe.second.op(i).has(eps)) {
                if(is_a<numeric>(fe.second.op(i)) && ex_to<numeric>(fe.second.op(i)).is_nonneg_integer()) {
                    rem *= pow(fe.first.op(i), fe.second.op(i));
                    continue;
                }
            }
            klst.append(fe.first.op(i));
            vlst.append(fe.second.op(i));
        }
        
        rem = mma_collect(rem, x(wild()));
        
        lst rem_lst;
        if(is_a<add>(rem)) {
            for(auto item : rem) rem_lst.append(item);
        } else {
            rem_lst.append(rem);
        }
        
        lst vlst2 = vlst;
        vlst2.append(1);
        for(auto item : rem_lst) {
            lst klst2 = klst;
            klst2.append(item);
            FunExp.push_back(make_pair(klst2, vlst2));
        }
    }
    
}

// Section 2.1 @ https://arxiv.org/pdf/1712.04441.pdf
// also refers to Feng/Thinking.pdf
void SD::Scalelesses(bool verb) {
    if(IsZero) return;
    if(Deltas.size()<1) return;
    if(verb) cout << now() << " - Scaleless: " << FunExp.size() << " :> " << flush;

    vector<ex> sl_res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, FunExp, [&](auto &funexp, auto rid) {
        symbol s;
        auto fun = funexp.first;
        auto exp = funexp.second;
        bool is0;
        for(auto delta : Deltas) {
            is0 = false;
            if(delta.nops()<2) continue;
            // to make sure the integrand is prejective
            if(delta.nops()>1) {
                lst sRepl;
                for(int j=0; j<delta.nops(); j++) {
                    sRepl.append(delta[j]==delta[j]*s);
                }
                
                bool is_s = true;
                ex n_s = 0;
                for(int j=0; j<fun.nops(); j++) {
                    auto tmp = fun.op(j).subs(sRepl).expand();
                    if(tmp.degree(s)!=tmp.ldegree(s)) {
                        is_s = false;
                        break;
                    }
                    n_s += tmp.degree(s) * exp.op(j);
                }
                if(!is_s || !normal(n_s+delta.nops()).is_zero()) continue;
            }
            
            for(long long i=1; i<ex_to<numeric>(GiNaC::pow(2,delta.nops())).to_long()-1; i++) {
                lst sRepl;
                auto ci = i;
                for(int j=0; (j<delta.nops() && ci>0); j++) {
                    if((ci % 2)==1) sRepl.append(delta[j]==delta[j]*s);
                    ci = ci / 2;
                }
                
                bool is_s = true;
                ex n_s = 0;
                for(int j=0; j<fun.nops(); j++) {
                    if(is_a<numeric>(exp.op(j)) && ex_to<numeric>(exp.op(j)).is_nonneg_integer()) continue;
                    auto tmp = mma_collect(fun.op(j).subs(sRepl),s);
                    if(tmp.degree(s)!=tmp.ldegree(s)) {
                        is_s = false;
                        break;
                    }
                    n_s += tmp.degree(s) * exp.op(j);
                }
                if(!is_s) continue;
                if(normal(n_s).has(ep)) {
                    is0 = true;
                    break;
                }
            }
            if(is0) break;
        }
        if(!is0) return lst{fun, exp};
        else return lst{};
    }, "SL", 0, true);
    
    FunExp.clear();
    FunExp.shrink_to_fit();
    for(auto item : sl_res) {
        lst kv = ex_to<lst>(item);
        if(kv.nops()<1) continue;
        FunExp.push_back(make_pair(ex_to<lst>(kv.op(0)), ex_to<lst>(kv.op(1))));
    }
    
    if(verb) cout << FunExp.size() << endl;
    if(FunExp.size()<1) IsZero = true;
}

void SD::RemoveDeltas() {
    if(IsZero) return;
    if(Deltas.size()<1) return;
    
    vector<pair<lst,lst>> funexp = FunExp;
    for(auto xs : Deltas) {
        vector<pair<lst,lst>> tmp;
        for(int i=0; i<xs.nops(); i++) {
            auto xj = xs.op(i);
            auto jInv = lstHelper::sum(xs).subs(xj==1);
            
            exmap repl;
            for(int j=0; j<xs.nops(); j++) {
                auto xxj = xs.op(j);
                if(xxj != xj) repl[xxj] = xj*xxj;
            }
            
            for(auto fe : funexp) {
                lst funs;
                lst exps = fe.second;
                ex expns = 0;
                for(int j=0; j<fe.first.nops(); j++) {
                    auto fun = fe.first.op(j);
                    fun = fun.subs(repl).normal();
                    if(!fun.is_polynomial(xj)) {
                        cerr << "xj: " << xj << endl;
                        cerr << "fun: " << fun << endl;
                        cerr << "fun is NOT polynormial of xj." << endl;
                        assert(false);
                    }
                    auto expn = mma_collect(fun, xj).degree(xj);
                    fun = pow(xj, -expn) * fun;
                    fun = normal(fun.subs(xj==1/xj));
                    fun = fun.subs(xj==jInv);
                    funs.append(fun);
                    expns += expn * exps.op(j);
                }
                
                funs.append(jInv);
                exps.append(ex(0)-xs.nops()-expns);
                tmp.push_back(make_pair(funs, exps));
            }
        }
        funexp = tmp;
    }
    FunExp = funexp;
    Deltas.clear();
    XReOrders();
    Normalizes();
}

// after SDPrepares, Integrands can be expanded in ep safely.
void SD::SDPrepares() {
    if(IsZero) return;
    if(FunExp.size()<1) {
        IsZero = true;
        return;
    }
    
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    //check vs
    for(auto &kv : FunExp) {
        ex ft = kv.first.op(1);
        if(ft.has(vs)) {
            ft = mma_collect(ft, vs);
            assert(ft.is_polynomial(vs) && (ft.degree(vs)-1)==0);
            ex expn = -kv.second.op(1);
            // (2*Pi*I) dropped out, since we will take residue later.
            kv.first.let_op(0) = kv.first.op(0) * tgamma(expn+vz)*tgamma(-vz)/tgamma(expn)*pow(vs,vz);
            ex w1 = ft.coeff(vs);
            ex w2 = ft.subs(vs==0);
            if(!w2.is_zero()) {
                kv.first.let_op(1) = w2;
                kv.first.append(w1);
                kv.second.let_op(1) = kv.second.op(1)-vz;
                kv.second.append(vz);
            }
        }
    }
    
    Integrands.clear();
    Integrands.shrink_to_fit();
    
    if(CheckF1) {
        if(Verbose > 0) cout << now() << " - Bisection: " << FunExp.size() << " :> " << flush;
        vector<ex> funexps =
        GiNaC_Parallel(ParallelProcess, ParallelSymbols, FunExp, [&](auto &kv, auto rid) {
            lst para_res_lst;
            auto kvs = AutoF1(kv);
            for(auto item : kvs) {
                para_res_lst.append(lst{item.first, item.second});
            }
            return para_res_lst;
        }, "f1", 0, !debug);

        FunExp.clear();
        for(auto &item : funexps) {
            for(auto &it : ex_to<lst>(item)) {
                FunExp.push_back(make_pair(ex_to<lst>(it.op(0)), ex_to<lst>(it.op(1))));
            }
        }
        if(Verbose > 0) cout << FunExp.size() << endl;
    }
    
    auto kvs = FunExp;
    FunExp.clear();
    FunExp.shrink_to_fit();
    for(auto &kv : kvs) {
        bool to_add = true;
        for(auto item : kv.first) {
            if(item.is_zero()) {
                to_add = false;
                break;
            }
        }
        if(to_add) FunExp.push_back(kv);
    }
    if(FunExp.size()<1) {
        IsZero = true;
        return;
    }
    
    if(Verbose > 0) cout << now() << " - SDPrepares: ..." << endl << flush;
    
    vector<ex> sd_res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, FunExp, [&](auto &kv, auto rid) {
        // return a lst, element pattern: { {{x1,n1}, {x2,n2}, ...}, {{e1, n1},{e2,n2}, ...} }.
        lst para_res_lst;
        auto xns_pns = DS(kv);
        for(auto const &item : xns_pns) {

            //TODO: handle other poles, e.g., PolyGamma/Gamma from input
            // take z-poles
            if(item.has(vz)) {
                lst zpols;
                // poles from Gamma(-z)
                for(int vn=0; vn<=sN; vn++) {
                    zpols.append(vn);
                }
                // poles from xi^{c1*z+c0}
                for(auto xn : item.op(0)) {
                    assert(is_a<numeric>(xn.op(1).coeff(vz)));
                    if(xn.op(1).coeff(vz)<0) {
                        ex c0 = xn.op(1).coeff(vz, 0);
                        ex c1 = xn.op(1).coeff(vz, 1);
                        int pxn = -1;
                        while(true) {
                            ex zp = (pxn-c0)/c1;
                            ex zpn = zp.subs(lst{eps==0,ep==0});
                            assert(is_a<numeric>(zpn));
                            if(zpn>sN) break;
                            zpols.append(zp);
                            pxn -= 1;
                        }
                    }
                }
                zpols.sort();
                zpols.unique();
                
                symbol ss;
                for(auto zp : zpols) {
                    auto xn = item.op(0);
                    auto pn = item.op(1);
                    for(int i=0; i<xn.nops(); i++) {
                        xn.let_op(i).let_op(0) = xn.op(i).op(0).subs(vz==ss+zp).subs(ss==vz);
                        xn.let_op(i).let_op(1) = xn.op(i).op(1).subs(vz==ss+zp).subs(ss==vz);
                    }
                    for(int i=0; i<pn.nops(); i++) {
                        pn.let_op(i).let_op(0) = pn.op(i).op(0).subs(vz==ss+zp).subs(ss==vz);
                        pn.let_op(i).let_op(1) = pn.op(i).op(1).subs(vz==ss+zp).subs(ss==vz);
                    }
                    para_res_lst.append(lst{xn, pn});
                }
            } else {
                para_res_lst.append(item);
            }
        }
        return para_res_lst;
    }, "SD", Verbose, true);
    
    ex min_expn = 1;
    vector<ex> ibp_in_vec;
    for(auto &item : sd_res) {
        for(auto &it : ex_to<lst>(item)) {
            ex expn = 0;
            for(auto xn : it.op(0)) {
                ex nxn = xn.op(1).subs(lst{ep==0, eps==0, vz==0});
                if(nxn<-1) expn += nxn+1;
            }
            if(expn < min_expn) min_expn = expn;
            
            int sim_max;
            if((ex(0)-expn)>=10) sim_max = 1;
            else if((ex(0)-expn)>=8) sim_max = 3;
            else if((ex(0)-expn)>=6) sim_max = 5;
            else if((ex(0)-expn)>=4) sim_max = 10;
            else if((ex(0)-expn)>=2) sim_max = 50;
            else sim_max = 100;
            
            //sim_max = 0; //disable sim
                        
            lst xns = ex_to<lst>(it.op(0));
            lst pns;
            ex tmp = 1;
            for(int i=0; i<it.op(1).nops(); i++) {
                lst pn = ex_to<lst>(it.op(1).op(i));
                bool sim = pn.op(0).expand().nops()<=sim_max;
                bool nni = is_a<numeric>(pn.op(1)) && ex_to<numeric>(pn.op(1)).is_nonneg_integer();
                if(i>1 && (sim||nni)) {
                    tmp *= pow(pn.op(0), pn.op(1));
                } else if(i<2 || pn.op(0)!=1) {
                    pns.append(pn);
                }
            }
            if(tmp!=1) pns.append(lst{tmp, 1});
            ibp_in_vec.push_back(lst{xns, pns});
        }
    }
    if(Verbose > 1) cout << WHITE << "  \\--Maximum x^n: " << ex(0)-min_expn << "+1" << RESET << endl << flush;

    int pn = 0;
    vector<ex> ibp_res_vec;
    while(ibp_in_vec.size()>0) {
        pn++;
        ostringstream spn;
        spn << "IBP-" << pn;
        vector<ex> ibp_res =
        GiNaC_Parallel(ParallelProcess, ParallelSymbols, ibp_in_vec, [&](auto &xns_pns, auto rid) {
            // return lst
            // if(the lst length is 1): { 1 element for pole reached }
            // else: { elements for pole NOT reached, need to loop again }.
            // element pattern still as { {{x1,n1}, {x2,n2}, ...}, {{e1, n1},{e2,n2}, ...} }
            
            auto xns = xns_pns.op(0);
            auto pns = xns_pns.op(1);
            
            exset fts;
            pns.op(0).find(FTX(wild(1),wild(2)), fts);
            bool noFT = (fts.size()==1) && ( (*(fts.begin())).op(0) == 1 );
            
            ex pole_requested = -1;
            if(noFT || PoleRequested > -1) pole_requested = PoleRequested;
            
            for(int n=0; n<xns.nops(); n++) {
                ex xn = xns.op(n);
                auto expn = xn.op(1).subs(lst{eps==0,ep==0,vz==0}).normal();
                if(!is_a<numeric>(expn)) {
                    cout << RED << "expn NOT numeric: " << expn << RESET << endl;
                    assert(false);
                }

                if(ex_to<numeric>(expn) < pole_requested) {
                    auto xx = xn.op(0);
                    pns.let_op(0).let_op(0) = pns.op(0).op(0) / (xn.op(1)+1);
                    
                    lst xns2;
                    for(int i=0; i<xns.nops(); i++) {
                        if(i!=n) xns2.append(xns.op(i));
                    }
                    lst pns2 = ex_to<lst>(pns);
                    for(int i=0; i<pns.nops(); i++) {
                        pns2.let_op(i).let_op(0) = pns2.op(i).op(0).subs(xx==1);
                    }
                    
                    lst xns_pns_lst;
                    xns_pns_lst.append(lst{xns2, pns2});
                    
                    lst xns3 = ex_to<lst>(xns);
                    xns3.let_op(n).let_op(1) = xn.op(1)+1;
                    
                    for(int i=0; i<pns.nops(); i++) {
                        ex tmp = ex(0)-pns.op(i).op(1)*mma_diff(pns.op(i).op(0),xx,1,false);
                        if(tmp.is_zero()) continue;
                        
                        auto xs = get_x_from(tmp);
                        bool tz = false;
                        for(auto xi : xs) {
                            if(tmp.subs(xi==0).is_zero()) {
                                tz = true;
                                break;
                            }
                        }
                        // need collect_common_factors
                        if(tz) tmp = collect_common_factors(tmp);
                        if(tmp.is_zero()) continue;
                        
                        if(tz && is_a<mul>(tmp)) {
                            ex rem = 1;
                            for(auto ii : tmp) {
                                if(ii.match(pow(x(wild()), wild(2)))) {
                                    bool t = true;
                                    for(int ij=0; ij<xns3.nops(); ij++) {
                                        if(xns3.op(ij).op(0)==ii.op(0)) {
                                            xns3.let_op(ij).let_op(1) += ii.op(1);
                                            t = false;
                                            break;
                                        }
                                    }
                                    if(t) xns3.append(lst{ii.op(0), ii.op(1)});
                                } else if(ii.match(x(wild()))) {
                                    bool t = true;
                                    for(int ij=0; ij<xns3.nops(); ij++) {
                                        if(xns3.op(ij).op(0)==ii) {
                                            xns3.let_op(ij).let_op(1) += 1;
                                            t = false;
                                            break;
                                        }
                                    }
                                    if(t) xns3.append(lst{ii, 1});
                                } else {
                                    rem *= ii;
                                }
                            }
                            tmp = rem;
                        }
                        
                        lst pns3 = ex_to<lst>(pns);
                        if(pns.op(i).op(1)==1) {
                            pns3.let_op(i).let_op(0) = tmp;
                        } else {
                            pns3.let_op(i).let_op(1) = pns.op(i).op(1)-1;
                            int nn = pns.nops();
                            if(!(pns.op(nn-1).op(1)-1).is_zero()) pns3.append(lst{ tmp, 1 });
                            else pns3.let_op(nn-1).let_op(0) = pns.op(nn-1).op(0) * tmp;
                        }
                        xns_pns_lst.append(lst{xns3, pns3});
                    }
                    
                    return xns_pns_lst;
                }
            }
            
            return lst { xns_pns };

        }, spn.str().c_str(), Verbose, true);
    
        ibp_in_vec.clear();
        ibp_in_vec.shrink_to_fit();
        for(auto &item : ibp_res) {
            if(item.nops()>1) {
                for(auto &it : ex_to<lst>(item)) ibp_in_vec.push_back(it);
            } else {
                ex expr = 1;
                for(auto pn : item.op(0).op(1)) expr *= pow(pn.op(0), pn.op(1));
                ibp_res_vec.push_back(lst{ item.op(0).op(0), expr });
            }
        }
    }
    
    vector<ex> res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, ibp_res_vec, [&](auto &xns_expr, auto rid) {

        // return single element in which ep/eps can be expanded safely.
        lst para_res_lst;
        auto xns = xns_expr.op(0);
        auto expr = xns_expr.op(1);
        lst exprs = { expr };
        symbol dx;
        for(auto xn : xns) {
            auto expn = xn.op(1).subs(lst{eps==0,ep==0,vz==0}).normal();
            assert(is_a<numeric>(expn));
            
            lst exprs2;
            for(auto it : exprs) {
                ex rem = pow(xn.op(0), xn.op(1)) * it;
                if(ex_to<numeric>(expn)<=-1) {
                    ex dit = it;
                    ex dit0 = dit.subs(xn.op(0)==0);
                    ex ifact = 1;
                    rem -= pow(xn.op(0), xn.op(1)) * dit0 / ifact;
                    exprs2.append(dit0/(xn.op(1)+1)/ifact);
                    for(int i=1; i+expn<0; i++) {
                        dit = mma_diff(dit, xn.op(0));
                        dit0 = dit.subs(xn.op(0)==0);
                        ifact *= i;
                        rem -= pow(xn.op(0), xn.op(1)+i) * dit0 / ifact;
                        exprs2.append(dit0/(xn.op(1)+i+1)/ifact);
                    }
                }
                exprs2.append(rem);
            }
            exprs = exprs2;
        }

        for(auto const &it : exprs) {
            if(!it.is_zero()) para_res_lst.append(it);
        }
     
        for(int i=0; i<para_res_lst.nops(); i++) {
            auto xs = get_x_from(para_res_lst.op(i));
            
            lst x2y;
            for(int i=0; i<xs.size(); i++) {
                x2y.append(xs[i]==y(i));
            }
            
            para_res_lst.let_op(i) = para_res_lst.op(i).subs(x2y).subs(y(wild())==x(wild()));
        }

        //deleted from GiNaC 1.7.7
        //if(para_res_lst.nops()<1) para_res_lst.append(0);
        return para_res_lst;
    }, "Taylor", Verbose, true);
    
    // Take z-residues
    bool zResides = false;
    vector<ex> ints;
    for(auto &item : res) {
        for(auto it : ex_to<lst>(item)) {
            if(!it.is_zero()) ints.push_back(it);
            if(!zResides && it.has(vz)) zResides = true;
        }
    }
    
    if(zResides) {
        Integrands =
        GiNaC_Parallel(ParallelProcess, ParallelSymbols, ints, [&](auto &item, auto rid) {
            ex it = item;
            if(it.has(vz)) {
                exset cts;
                it.find(CT(wild()),cts);
                lst repl;
                for(auto ct : cts) {
                    ex cc = 1;
                    ex ll = 1;
                    lst cls;
                    if(is_a<mul>(ct.op(0))) {
                        for(auto ii : ct.op(0)) cls.append(ii);
                    } else {
                        cls.append(ct.op(0));
                    }
                    for(auto cl : cls) {
                        if(cl.has(vz)) cc *= cl;
                        else ll *= cl;
                    }
                    if(cc!=1) repl.append(ct==cc*CT(ll));
                }
                it = it.subs(repl);
                it = mma_series(it,vz,-1);
                it = ex(0)-it.coeff(vz, -1);
            }
            return it;
        }, "zResidue", Verbose, true);
    } else {
        Integrands = ints;
    }

}

void SD::EpsEpExpands() {
    if(IsZero) return;
    if(Integrands.size()<1) {
        IsZero = true;
        return;
    }
    
    if(Verbose > 0) cout << now() << " - EpsEpExpands ..." << endl << flush;
    
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    vector<ex> res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, Integrands, [&](auto &item, auto rid) {
        // return { {two elements}, {two elements}, ...},
        // 1st: x-independent coefficient, expanded in ep/eps
        // 2nd: x-integrand
        if(item.is_zero()) return lst{ lst{0, 0} };
        exset cts;
        item.find(CT(wild()), cts);
        if(cts.size() != 1) {
            cerr << "item: " << item << endl;
            cerr << "CT size is NOT 1: " << cts << endl;
            assert(false);
        }
        ex ct = (*(cts.begin())).subs(CT(wild())==wild()).subs(iEpsilon==0);
        auto it = item.subs(CT(wild())==1);
        it = mma_collect(it, vs, true);
        lst its;
        if(is_a<add>(it)) {
            for(auto ii : it) its.append(ii);
        } else {
            its.append(it);
        }

        lst para_res_lst;
        for(int i=0; i<its.nops();i++) {
            auto tmp = its.op(i);
            auto vc = tmp.subs(CCF(wild())==1);
            tmp = tmp / vc;
            tmp = tmp.subs(CCF(wild())==wild());
            //if(use_CCF) tmp = collect_common_factors(tmp);
            if(!tmp.has(eps) && !ct.has(eps)) {
                int ctN = epRank(ct);
                tmp = mma_series(tmp, ep, epN-ctN);
                for(int di=tmp.ldegree(ep); (di<=tmp.degree(ep) && di<=epN-ctN); di++) {
                    auto intg = tmp.coeff(ep, di);
                    assert(!intg.has(ep));
                    auto pref = mma_series(ct, ep, epN-di);
                    //if(use_CCF) intg = collect_common_factors(intg);
                    para_res_lst.append(lst{pref * pow(ep, di) * vc, intg});
                }
            } else {
                auto sct = ct;
                int sctN = epsRank(sct);
                ex stmp = mma_series(tmp, eps, epsN-sctN);
                for(int sdi=stmp.ldegree(eps); (sdi<=stmp.degree(eps) && sdi<=epsN-sctN); sdi++) {
                    tmp = stmp.coeff(eps, sdi);
                    //if(use_CCF) tmp = collect_common_factors(tmp);
                    assert(!tmp.has(eps));
                    ct = mma_series(sct, eps, epsN-sdi);
                    int ctN = epRank(ct);
                    tmp = mma_series(tmp, ep, epN-ctN);
                    for(int di=tmp.ldegree(ep); (di<=tmp.degree(ep) && di<=epN-ctN); di++) {
                        auto intg = tmp.coeff(ep, di);
                        assert(!intg.has(ep));
                        auto pref = mma_series(ct, ep, epN-di);
                        //if(use_CCF) intg = collect_common_factors(intg);
                        para_res_lst.append(lst{pref * pow(eps, sdi) * pow(ep, di) * vc, intg});
                    }
                }
            }
        }
        
        //deleted from GiNaC 1.7.7
        //if(para_res_lst.nops()<1) para_res_lst.append(lst{0,0});
        return para_res_lst;

    }, "EpsEp", Verbose, !debug);
    
    if(Verbose > 1) cout << "  \\--Collecting: ";
    map<ex, ex, ex_is_less> int_pref;
    long long ncollect = 0;
    for(auto &item : res) {
        ncollect += item.nops();
        for(auto &kv : ex_to<lst>(item)) {
            int_pref[kv.op(1)] += kv.op(0);
        }
    }
    
    if(Verbose > 1) cout << ncollect << " :> " << flush;
    expResult.clear();
    expResult.shrink_to_fit();
    for(auto kv : int_pref) {
        if(kv.second.normal().is_zero()) continue;
        expResult.push_back(make_pair(kv.second, kv.first));
    }
    if(Verbose > 1) cout << expResult.size() << endl;
}

void SD::CompileMatDet() {
    auto pid = getpid();
    std::ofstream ofs;
    ostringstream cppfn, cmd;
    cppfn << pid << "/MatDet.cpp";
    ofs.open(cppfn.str(), ios::out);
    if (!ofs) throw runtime_error("failed to open *.cpp file!");
/*----------------------------------------------*/
ofs << R"EOF(
#include <math.h>
#include <complex>
extern "C" {
#include <quadmath.h>
}
#define Pi 3.1415926535897932384626433832795028841971693993751L
#define Euler 0.57721566490153286060651209008240243104215933593992L

using namespace std;
typedef __float128 qREAL;
typedef __complex128 qCOMPLEX;
typedef long double dREAL;
typedef complex<long double> dCOMPLEX;

dREAL expt(dREAL a, dREAL b) { return pow(a,b); }
dCOMPLEX expt(dCOMPLEX a, dREAL b) { return pow(a,b); }
dREAL recip(dREAL a) { return 1.L/a; }
dCOMPLEX recip(dCOMPLEX a) { return 1.L/a; }

qREAL expt(qREAL a, qREAL b) { return powq(a,b); }
qCOMPLEX expt(qCOMPLEX a, qREAL b) { return cpowq(a,b); }
qREAL recip(qREAL a) { return 1.Q/a; }
qCOMPLEX recip(qCOMPLEX a) { return 1.Q/a; }

qREAL pow(qREAL x, qREAL y) { return powq(x, y); }
qREAL log(qREAL x) { return logq(x); }
qCOMPLEX pow(qCOMPLEX x, qREAL y) { return cpowq(x, y); }
qCOMPLEX log(qCOMPLEX x) { return clogq(x); }

dCOMPLEX MatDetL(dCOMPLEX mat[], int n) {
    bool is_zero = false;
    int s=1;
    for(int i=0; i<n-1; i++) {
        if(fabs(mat[i*n+i])<1.0E-15) {
            bool is_zero = true;
            for(int j=i+1; j<n; j++) {
                if(fabs(mat[i*n+j])>1.0E-15) {
                    for(int k=0; k<n; k++) {
                        auto tmp = mat[k*n+j];
                        mat[k*n+j] = mat[k*n+i];
                        mat[k*n+i] = tmp;
                    }
                    is_zero = false;
                    s=-s;
                    break;
                }
            }
            if(is_zero) return 0;
        }
        for(int k=i+1; k<n; k++) {
            auto m = mat[k*n+i]/mat[i*n+i];
            for(int j=0; j<n; j++) mat[k*n+j] = mat[k*n+j] - m*mat[i*n+j];
        }
    }
    dCOMPLEX ret = s;
    for(int k=0; k<n; k++) ret *= mat[k*n+k];
    return ret;
}

#undef Pi
#undef Euler
#define Pi 3.1415926535897932384626433832795028841971693993751Q
#define Euler 0.57721566490153286060651209008240243104215933593992Q

qCOMPLEX MatDetQ(qCOMPLEX mat[], int n) {
    bool is_zero = false;
    int s=1;
    for(int i=0; i<n-1; i++) {
        if(cabsq(mat[i*n+i])<1.0E-15) {
            bool is_zero = true;
            for(int j=i+1; j<n; j++) {
                if(cabsq(mat[i*n+j])>1.0E-15) {
                    for(int k=0; k<n; k++) {
                        auto tmp = mat[k*n+j];
                        mat[k*n+j] = mat[k*n+i];
                        mat[k*n+i] = tmp;
                    }
                    is_zero = false;
                    s=-s;
                    break;
                }
            }
            if(is_zero) return 0;
        }
        for(int k=i+1; k<n; k++) {
            auto m = mat[k*n+i]/mat[i*n+i];
            for(int j=0; j<n; j++) mat[k*n+j] = mat[k*n+j] - m*mat[i*n+j];
        }
    }
    qCOMPLEX ret = s;
    for(int k=0; k<n; k++) ret *= mat[k*n+k];
    return ret;
}
)EOF" << endl;
/*----------------------------------------------*/
    ofs.close();
    cmd.clear();
    cmd.str("");
    cmd << "g++ -fPIC " << CFLAGS << " -c -o " << pid << "/MatDet.o " << pid << "/MatDet.cpp";
    system(cmd.str().c_str());
}

void SD::CIPrepares(const char *key) {
    if(IsZero) return;
    if(expResult.size()<1) {
        IsZero = true;
        return;
    }
    
    if(Verbose > 0) cout << now() << " - CIPrepares ..." << endl << flush;
    auto pid = getpid();
    
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    vector<ex> resf =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, expResult, [&](auto &kv, auto rid) {
        // return lst{ kv.first, kv.second, ft};
        auto expr = kv.second;
        auto xs = get_xy_from(expr);
        if(xs.size()<1) {
            return lst{kv.first, kv.second, 1};
        }

        exset ftxset;
        expr.find(FTX(wild(1),wild(2)), ftxset);
        ex ft;
        int ftxsize = -1;
        for(auto item : ftxset) {
            auto xys = get_xy_from(item.op(0));
            if((int)xys.size() > ftxsize) {
                ft = item.op(0);
                ftxsize = xys.size();
            }
        }
        
        bool need_contour_deformation = false;
        if(ft.has(x(wild())) && !ft.has(PL(wild()))) {
            auto tmp = ft.subs(nReplacements).expand();
            if(is_a<add>(tmp)) {
                for(auto item : tmp) {
                    assert(is_a<numeric>(item.subs(x(wild())==1)));
                    if(item.subs(x(wild())==1) < 0) {
                        need_contour_deformation = true;
                        break;
                    }
                }
            } else {
                assert(is_a<numeric>(tmp.subs(x(wild())==1)));
                if(tmp.subs(x(wild())==1) < 0) need_contour_deformation = true;
            }
            if(!need_contour_deformation) ft = 1; //note the difference with SDPrepare
        } else if(!ft.has(x(wild()))){
            ft = 1;
        } else {
            double tmin = FindMinimum(ft, true);
            if(tmin > 1E-5) {
                ft = 1;
            } 
        }
        
        ft = collect_common_factors(ft);
        return lst{ kv.first, kv.second, ft};
        
    }, "CI-F", Verbose, false);
    

//============================================================================================================
    lst fts;
    for(auto item : resf) {
        if(item.op(2).has(x(wild()))) {
            fts.append(item.op(2));
        }
    }
    fts.sort();
    fts.unique();
    
    vector<pair<ex,int>> ftnvec;
    map<ex,int,ex_is_less> ftnmap;
    int ft_n = 1;
    FT_N_NX.remove_all();
    
    //deleted from GiNaC 1.7.7
    //FT_N_NX.append(lst{0, 0});
    
    for(auto item : fts) {
        ftnvec.push_back(make_pair(item, ft_n));
        ftnmap[item] = ft_n;
        FT_N_NX.append(lst{ft_n, get_xy_from(item).size()});
        ft_n++;
    }
    //ftnvec item: lst { ft, ft-id }
    
    vector<lst> res_vec;
    map<ex, ex, ex_is_less> cf_int;
    
    for(auto &item : resf) {
        auto ii = ex_to<lst>(item);
        if(ii.op(2)==1) {
            ii.append(-1);
        } else {
            int ft_n = ftnmap[item.op(2)];
            if(ft_n==0) {
                cerr << item.op(2) << endl;
                assert(false);
            }
            ii.append(ft_n);
        }
        if(ii.op(0).is_zero() || ii.op(1).subs(FTX(wild(1),wild(2))==1).is_zero()) continue;
        if(!use_IBF) res_vec.push_back(ii);
        else {
            ex key = ii;
            key.let_op(1) = 1;
            cf_int[key] += ii.op(1);
        }
    }
    
    if(use_IBF) {
        for(auto kv : cf_int) {
            lst ii = ex_to<lst>(kv.first);
            ii.let_op(1) = kv.second;
            res_vec.push_back(ii);
        }
    }
    
    //res_vec item: lst { coeff, integrand, ft, ft-id }
    
    if(res_vec.size()<1) {
        IsZero = true;
        return;
    }
//============================================================================================================

    // Prepare FT-lambda
    if(use_cpp)
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, ftnvec, [&](auto &kv, auto rid) {
        // return nothing
        ex ft = kv.first;
        ex ft_n = kv.second;
        auto fxs = get_xy_from(ft);
        lst las;
        
        // exp-fn configure, Exp[-Power[ftn,2*f2n]]//N
        char * ri = "0.0";
        int ftn = 10, f2n = 1;
        bool use_exp = true;
        ex ft_max = ex(0)-numeric(FindMinimum(-abs(ft)));
        ft_max = ft_max/ex(ftn);
        
        auto pls = get_pl_from(ft);
        int npls = pls.size()>0 ? ex_to<numeric>(pls[pls.size()-1].subs(lst{PL(wild())==wild()})).to_int() : -1;
        lst plRepl;
        for(int i=0; i<npls+1; i++) {
            ostringstream pl;
            pl << "pl[" << i << "]";
            plRepl.append(PL(i) == symbol(pl.str()));
        }
        
        ex DFs[fxs.size()], DDFs[fxs.size()*fxs.size()];
        for(int i=0; i<fxs.size(); i++) {
            auto df = mma_diff(ft, fxs[i]);
            DFs[i] = collect_common_factors(df);
            ostringstream ilaos;
            ilaos << "ilas[" << i << "]";
            symbol ila(ilaos.str());
            for(int j=0; j<fxs.size(); j++) {
                auto ddf = mma_diff(DFs[i], fxs[j]);
                DDFs[fxs.size()*i+j] = collect_common_factors(ddf);
            }
        }

        ostringstream cppfn, sofn;
        cppfn << pid << "/F" << ft_n << ".cpp";
        sofn << pid << "/F" << ft_n << ".o";
        std::ofstream ofs;
        ofs.open(cppfn.str(), ios::out);
        if (!ofs) throw runtime_error("failed to open *.cpp file!");
        
        lst cxRepl, czRepl;
        for (int i=0; i<fxs.size(); i++) {
            ostringstream sx, sz;
            sx << "x[" << i << "]";
            cxRepl.append(fxs[i] == symbol(sx.str()));
            sz << "z[" << i << "]";
            czRepl.append(fxs[i] == symbol(sz.str()));
        }

/*----------------------------------------------*/
ofs << R"EOF(
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <iostream>
extern "C" {
#include <quadmath.h>
}

using namespace std;

#define Pi 3.1415926535897932384626433832795028841971693993751L
#define Euler 0.57721566490153286060651209008240243104215933593992L

typedef __float128 qREAL;
typedef __complex128 qCOMPLEX;
typedef long double dREAL;
typedef complex<long double> dCOMPLEX;

dREAL expt(dREAL a, dREAL b);
dCOMPLEX expt(dCOMPLEX a, dREAL b);
dREAL recip(dREAL a);
dCOMPLEX recip(dCOMPLEX a);

qREAL expt(qREAL a, qREAL b);
qCOMPLEX expt(qCOMPLEX a, qREAL b);
qREAL recip(qREAL a);
qCOMPLEX recip(qCOMPLEX a);

qREAL pow(qREAL x, qREAL y);
qREAL log(qREAL x);
qCOMPLEX pow(qCOMPLEX x, qREAL y);
qCOMPLEX log(qCOMPLEX x);

)EOF" << endl;
/*----------------------------------------------*/
        auto cppL = CppFormat(ofs, "L");
        auto cppQ = CppFormat(ofs, "Q");
        ft = Evalf(ft);
        
        // FL_fid
        ofs << "dREAL FL_" << ft_n << "(const dREAL* x, const dREAL *pl) {" << endl;
        ofs << "dREAL yy = " << endl;
        ft.subs(plRepl).subs(cxRepl).print(cppL);
        ofs << ";" << endl;
        ofs << "return yy;" << endl; // find max image part, check with 0
        ofs << "}" << endl;
        ofs << endl;
        
        // FQ_fid
        ofs << "qREAL FQ_" << ft_n << "(const qREAL* x, const qREAL *pl) {" << endl;
        ofs << "qREAL yy = " << endl;
        ft.subs(plRepl).subs(cxRepl).print(cppQ);
        ofs << ";" << endl;
        ofs << "return yy;" << endl; // find max image part, check with 0
        ofs << "}" << endl;
        ofs << endl;
        
        // D's FL_fid
        ofs << "dREAL FL_" << ft_n << "(const int i, const dREAL* x, const dREAL *pl) {" << endl;
        for(int i=0; i<fxs.size(); i++) {
            ofs << "if("<<i<<"==i) return ";
            DFs[i].subs(plRepl).subs(cxRepl).print(cppL);
            ofs << ";" << endl;
        }
        ofs << "return 0;" << endl;
        ofs << "}" << endl;
        ofs << endl;
        
        // D's FQ_fid
        ofs << "qREAL FQ_" << ft_n << "(const int i, const qREAL* x, const qREAL *pl) {" << endl;
        for(int i=0; i<fxs.size(); i++) {
            ofs << "if("<<i<<"==i) return ";
            DFs[i].subs(plRepl).subs(cxRepl).print(cppQ);
            ofs << ";" << endl;
        }
        ofs << "return 0;" << endl;
        ofs << "}" << endl;
        ofs << endl;
        
        // DD's FL_fid
        ofs << "dREAL FL_" << ft_n << "(const int i, const int j, const dREAL* x, const dREAL *pl) {" << endl;
        for(int i=0; i<fxs.size(); i++) {
        for(int j=0; j<fxs.size(); j++) {
            ofs << "if("<<i<<"==i && "<<j<<"==j) return ";
            DDFs[i*fxs.size()+j].subs(plRepl).subs(cxRepl).print(cppL);
            ofs << ";" << endl;
        }}
        ofs << "return 0;" << endl;
        ofs << "}" << endl;
        ofs << endl;
        
        // DD's FQ_fid
        ofs << "qREAL FQ_" << ft_n << "(const int i, const int j, const qREAL* x, const qREAL *pl) {" << endl;
        for(int i=0; i<fxs.size(); i++) {
        for(int j=0; j<fxs.size(); j++) {
            ofs << "if("<<i<<"==i && "<<j<<"==j) return ";
            DDFs[i*fxs.size()+j].subs(plRepl).subs(cxRepl).print(cppQ);
            ofs << ";" << endl;
        }}
        ofs << "return 0;" << endl;
        ofs << "}" << endl;
        ofs << endl;
        
        // X2ZL_fid
        ofs << "void X2ZL_" << ft_n << "(const dREAL* x, dCOMPLEX* z, dCOMPLEX* r, dREAL* dff, const dREAL* pl, const dREAL* las) {" << endl;
        ofs << "int nfxs="<<fxs.size()<<", f2n="<<f2n<<";" << endl;
        ofs << "dREAL fmax = "; ft_max.print(cppL); ofs << ";" << endl;
        ofs << "dCOMPLEX ilas[nfxs];" << endl;
        ofs << "for(int i=0; i<nfxs; i++) ilas[i] = complex<long double>("<<ri<<"L*las[i], las[i]);" << endl;
        ofs << "dff[nfxs] = FL_"<<ft_n<<"(x,pl);" << endl;
        if(use_exp) ofs << "dREAL expf2n=exp(-pow(dff[nfxs]/fmax,2*f2n));" << endl;
        else ofs << "dREAL expf2n=1.L;" << endl;
        ofs << "for(int i=0; i<nfxs; i++) dff[i] = FL_"<<ft_n<<"(i,x,pl);" << endl;
        ofs << "for(int i=0; i<nfxs; i++) r[i] = dff[i]*ilas[i]*expf2n;" << endl;
        ofs << "for(int i=0; i<nfxs; i++) z[i] = x[i]-x[i]*(1-x[i])*r[i];" << endl;
        ofs << "}" << endl;
        ofs << endl;
        
        // X2ZQ_fid
        ofs << "void X2ZQ_" << ft_n << "(const qREAL* x, qCOMPLEX* z, qCOMPLEX* r, qREAL* dff, const qREAL* pl, const qREAL* las) {" << endl;
        ofs << "int nfxs="<<fxs.size()<<", f2n="<<f2n<<";" << endl;
        ofs << "qREAL fmax = "; ft_max.print(cppQ); ofs << ";" << endl;
        ofs << "qCOMPLEX ilas[nfxs];" << endl;
        ofs << "for(int i=0; i<nfxs; i++) ilas[i] = las[i] * ("<<ri<<"Q+1.Qi);" << endl;
        ofs << "dff[nfxs] = FQ_"<<ft_n<<"(x,pl);" << endl;
        if(use_exp) ofs << "qREAL expf2n=expq(-powq(dff[nfxs]/fmax,2*f2n));" << endl;
        else ofs << "dREAL expf2n = 1.Q;" << endl;
        ofs << "for(int i=0; i<nfxs; i++) dff[i] = FQ_"<<ft_n<<"(i,x,pl);" << endl;
        ofs << "for(int i=0; i<nfxs; i++) r[i] = dff[i]*ilas[i]*expf2n;" << endl;
        ofs << "for(int i=0; i<nfxs; i++) z[i] = x[i]-x[i]*(1-x[i])*r[i];" << endl;
        ofs << "}" << endl;
        ofs << endl;
        
        // MatL_id
        ofs << "void MatL_"<<ft_n<<"(dCOMPLEX* mat, const dREAL* x, const dREAL* dff, const dREAL* pl, const dREAL* las) {" << endl;
        ofs << "int nfxs="<<fxs.size()<<", f2n="<<f2n<<";" << endl;
        ofs << "dREAL fmax = "; ft_max.print(cppL); ofs << ";" << endl;
        ofs << "dCOMPLEX ilas[nfxs];" << endl;
        ofs << "for(int i=0; i<nfxs; i++) ilas[i] = complex<long double>("<<ri<<"L*las[i], las[i]);" << endl;
        if(use_exp) ofs << "dREAL expf2n = exp(-pow(dff[nfxs]/fmax,2*f2n));" << endl;
        else ofs << "dREAL expf2n = 1.L;" << endl;
        ofs << "for(int i=0; i<nfxs; i++) {" << endl;
        ofs << "for(int j=0; j<nfxs; j++) {" << endl;
        ofs << "int ij = i*nfxs+j;" << endl;
        ofs << "if(i!=j) mat[ij] = 0;" << endl;
        ofs << "else mat[ij] = 1.L-(1.L-2.L*x[i])*dff[i]*ilas[i]*expf2n;" << endl;
        ofs << "mat[ij] = mat[ij]-x[i]*(1-x[i])*FL_"<<ft_n<<"(i,j,x,pl)*ilas[i]*expf2n;" << endl;
        if(use_exp) {
            ofs << "dREAL df2n = -2*f2n*dff[j]/fmax*pow(dff[nfxs]/fmax,2*f2n-1);" << endl;
            ofs << "mat[ij] = mat[ij]-x[i]*(1-x[i])*dff[i]*ilas[i]*expf2n*df2n;" << endl;
        }
        ofs << "}}" << endl;
        ofs << "}" << endl;
        ofs << endl;
        
        // MatQ_fid
        ofs << "void MatQ_"<<ft_n<<"(qCOMPLEX *mat, const qREAL* x, const qREAL* dff, const qREAL *pl, const qREAL *las) {" << endl;
        ofs << "int nfxs="<<fxs.size()<<", f2n="<<f2n<<";" << endl;
        ofs << "qREAL fmax = "; ft_max.print(cppQ); ofs << ";" << endl;
        ofs << "qCOMPLEX ilas[nfxs];" << endl;
        ofs << "for(int i=0; i<nfxs; i++) ilas[i] = las[i] * ("<<ri<<"Q+1.Qi);" << endl;
        if(use_exp) ofs << "qREAL expf2n = expq(-powq(dff[nfxs]/fmax,2*f2n));" << endl;
        else ofs << "qREAL expf2n = 1.Q;" << endl;
        ofs << "for(int i=0; i<nfxs; i++) {" << endl;
        ofs << "for(int j=0; j<nfxs; j++) {" << endl;
        ofs << "int ij = i*nfxs+j;" << endl;
        ofs << "if(i!=j) mat[ij] = 0;" << endl;
        ofs << "else mat[ij] = 1.Q-(1.Q-2.Q*x[i])*dff[i]*ilas[i]*expf2n;" << endl;
        ofs << "mat[ij] = mat[ij]-x[i]*(1-x[i])*FQ_"<<ft_n<<"(i,j,x,pl)*ilas[i]*expf2n;" << endl;
        if(use_exp) {
            ofs << "qREAL df2n = -2*f2n*dff[j]/fmax*powq(dff[nfxs]/fmax,2*f2n-1);" << endl;
            ofs << "mat[ij] = mat[ij]-x[i]*(1-x[i])*dff[i]*ilas[i]*expf2n*df2n;" << endl;
        }
        ofs << "}}" << endl;
        ofs << "}" << endl;
        ofs << endl;
        
        // for Minimization of F(z)-image, x[xn-1] is the lambda
        ofs << "extern \"C\" " << endl;
        ofs << "dREAL imgF_"<<ft_n<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las_in) {" << endl;
        ofs << "dREAL las[xn-1];" <<endl;
        ofs << "for(int i=0; i<xn-1; i++) las[i] = las_in[i]*x[xn-1];" <<endl;
        ofs << "dCOMPLEX z[xn], r[xn];" << endl;
        ofs << "dREAL dff[xn+1];" << endl;
        ofs << "X2ZL_"<<ft_n<<"(x,z,r,dff,pl,las);" << endl;
        ofs << "dCOMPLEX zf = ";
        ft.subs(plRepl).subs(czRepl).print(cppL);
        ofs << ";" << endl;
        ofs << "return -zf.imag()/x[xn-1];" << endl; // find max image part, check with 0
        ofs << "}" << endl;
        ofs << endl;
        
        // for Minimization of DF-i
        for(int i=0; i<fxs.size(); i++) {
            ofs << "extern \"C\" " << endl;
            ofs << "dREAL dirF_"<<ft_n<<"_"<<i<<"(const int xn, const dREAL* x, const dREAL *pl, const dREAL *las) {" << endl;
            ofs << "int f2n="<<f2n<<";" << endl;
            ofs << "dREAL fmax = "; ft_max.print(cppL); ofs << ";" << endl;
            if(use_exp) ofs << "dREAL expf2n = exp(-pow(FL_"<<ft_n<<"(x,pl)/fmax,2*f2n));" << endl;
            else ofs << "dREAL expf2n = 1.L;" << endl;
            ofs << "dREAL yy = FL_"<<ft_n<<"("<<i<<", x, pl)*expf2n;" << endl;
            ofs << "return -fabs(yy);" << endl;
            ofs << "}" << endl;
            ofs << endl;
        }
        
        ostringstream cmd;
        cmd << "g++ -fPIC -c " << CFLAGS << " -o " << sofn.str() << " " << cppfn.str();
        system(cmd.str().c_str());
        
        if(!file_exists(sofn.str().c_str())) {
            cmd.clear();
            cmd.str("");
            cmd << "cp " << sofn.str() << " . ";
            system(cmd.str().c_str());
        }
        
        return 0;
    
    }, "CI-C", Verbose, false);


//============================================================================================================

    // Prepare Integrand
    vector<ex> res;
    if(use_cpp) res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, res_vec, [&](auto &kvf, auto rid) {
        // return lst{ no-x-result, xn, x-indepent prefactor, ft_n }
        // or     lst{ id(SD(D|Q)_id in .so), xn, x-indepent prefactor, ft_n }
        
        auto expr = kvf.op(1);
        auto xs = get_xy_from(expr);
        auto ft_n = kvf.op(3);
        bool hasF = (ft_n>0);
        
        if(xs.size()<1) {
            return lst{
                expr.subs(FTX(wild(1),wild(2))==1).subs(iEpsilon==I*power(10,-50)),
                xs.size(), kvf.op(0), -1
            };
        }
        
        auto ft = kvf.op(2);
        auto fxs = get_xy_from(ft);
        
        exset ftxset;
        expr.find(FTX(wild(1),wild(2)), ftxset);
        lst ftxlst;
        for(auto it : ftxset) ftxlst.append(it);
        expr = mma_collect(expr, FTX(wild(1),wild(2)));
        vector<pair<ex,ex>> ft_expr;
        for(auto item : ftxlst) {
            ft_expr.push_back(make_pair(item.op(1), expr.coeff(item)));
        }
        
        lst cxRepl, czRepl, plRepl;
        for (int i=0; i<fxs.size(); i++) {
            ostringstream xs, zs;
            xs << "x[" << i << "]";
            zs << "z[" << i << "]";
            cxRepl.append(fxs[i] == symbol(xs.str()));
            czRepl.append(fxs[i] == symbol(zs.str()));
        }
        int count = fxs.size();
        for(auto xi : xs) {
            auto xii = xi.subs(czRepl);
            if(xii == xi) {
                ostringstream xs, zs;
                xs << "x[" << count << "]";
                cxRepl.append(xi == symbol(xs.str()));
                czRepl.append(xi == symbol(xs.str()));
                count++;
            }
        }
        assert(count==xs.size());
        auto pls = get_pl_from(expr);
        int npls = pls.size()>0 ? ex_to<numeric>(pls[pls.size()-1].subs(lst{PL(wild())==wild()})).to_int() : -1;
        for(int i=0; i<npls+1; i++) {
            ostringstream pl;
            pl << "pl[" << i << "]";
            plRepl.append(PL(i) == symbol(pl.str()));
        }
        
        ostringstream cppfn;
        cppfn << pid << "/" << rid << ".cpp";
        std::ofstream ofs;
        ofs.open(cppfn.str(), ios::out);
        if (!ofs) throw runtime_error("failed to open *.cpp file!");

/*----------------------------------------------*/
ofs << R"EOF(
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <iostream>
extern "C" {
#include <quadmath.h>
}

using namespace std;

typedef __float128 qREAL;
typedef __complex128 qCOMPLEX;
typedef long double dREAL;
typedef complex<long double> dCOMPLEX;

dCOMPLEX MatDetL(dCOMPLEX mat[], int n);
qCOMPLEX MatDetQ(qCOMPLEX mat[], int n);

dREAL expt(dREAL a, dREAL b);
dCOMPLEX expt(dCOMPLEX a, dREAL b);
dREAL recip(dREAL a);
dCOMPLEX recip(dCOMPLEX a);

qREAL expt(qREAL a, qREAL b);
qCOMPLEX expt(qCOMPLEX a, qREAL b);
qREAL recip(qREAL a);
qCOMPLEX recip(qCOMPLEX a);

qREAL pow(qREAL x, qREAL y);
qREAL log(qREAL x);
qCOMPLEX pow(qCOMPLEX x, qREAL y);
qCOMPLEX log(qCOMPLEX x);

#define Pi 3.1415926535897932384626433832795028841971693993751L
#define Euler 0.57721566490153286060651209008240243104215933593992L
#define iEpsilon complex<long double>(0,1.E-50L)

)EOF" << endl;
/*----------------------------------------------*/

        if(hasF) {
            ofs << "dREAL FL_"<<ft_n<<"(const int, const dREAL*, const dREAL*);" << endl;
            ofs << "qREAL FQ_"<<ft_n<<"(const int, const qREAL*, const qREAL*);" << endl;
            ofs << "dREAL FL_"<<ft_n<<"(const int, const int, const dREAL*, const dREAL*);" << endl;
            ofs << "qREAL FQ_"<<ft_n<<"(const int, const int, const qREAL*, const qREAL*);" << endl;
            ofs << "void X2ZL_"<<ft_n<<"(const dREAL*, dCOMPLEX*, dCOMPLEX*, dREAL*, const dREAL*, const dREAL*);" << endl;
            ofs << "void X2ZQ_"<<ft_n<<"(const qREAL*, qCOMPLEX*, qCOMPLEX*, qREAL*, const qREAL*, const qREAL*);" << endl;
            ofs << "void MatL_"<<ft_n<<"(dCOMPLEX*, const dREAL*, const dREAL*, const dREAL*, const dREAL*);" << endl;
            ofs << "void MatQ_"<<ft_n<<"(qCOMPLEX*, const qREAL*, const qREAL*, const qREAL*, const qREAL*);" << endl;
            ofs << endl << endl;
        }

        // long double
        auto cppL =  CppFormat(ofs, "L");
        ofs << "extern \"C\" " << endl;
        ofs << "int SDD_"<<rid<<"(const unsigned int xn, const qREAL qx[], const unsigned int yn, qREAL y[], const qREAL qpl[], const qREAL qlas[]) {" << endl;
        ofs << "dREAL x[xn], x0[xn];" << endl;
        ofs << "for(int i=0; i<xn; i++) x[i] = qx[i];" << endl;
        ofs << "dREAL pl["<<(npls<0 ? 1 : npls+1)<<"];" << endl;
        ofs << "for(int i=0; i<"<<(npls+1)<<"; i++) pl[i] = qpl[i];" << endl;
        if(hasF) {
            ofs << "dCOMPLEX z[xn],r[xn];" << endl;
            ofs << "dREAL dff[xn+1];";
            ofs << "dCOMPLEX yy=0, ytmp, det;" << endl;
            ofs << "int ii, nfxs="<<fxs.size()<<";" << endl;
            ofs << "dREAL las[nfxs];" << endl;
            ofs << "for(int i=0; i<nfxs; i++) las[i] = qlas[i];" << endl;
            ofs << "dCOMPLEX mat[nfxs*nfxs];" << endl;
            for(auto &kv : ft_expr) {
                ofs << "{" << endl;
                if(SD::debug) {
                    ofs << "//debug-xs: " << kv.first << endl;
                    ofs << "//debug-int: " << kv.second << endl;
                }
                lst xs0;
                for(int ii=0; ii<fxs.size(); ii++) {
                    if(!kv.first.has(fxs[ii])) xs0.append(ii);
                }
                ofs << "for(int i=0; i<xn; i++) z[i] = x0[i] = x[i];" << endl;
                for(auto x0i : xs0) ofs << "x0["<<x0i<<"]=0;" << endl;
                ofs << "X2ZL_"<<ft_n<<"(x0,z,r,dff,pl,las);" << endl;
                ofs << "MatL_"<<ft_n<<"(mat,x0,dff,pl,las);" << endl;
                for(auto x0i : xs0) {
                    ofs << "ii = " << x0i << ";" << endl;
                    ofs << "z[ii] = x[ii]-x[ii]*(1-x[ii])*r[ii];" << endl;
                    ofs << "for(int j=0; j<nfxs;j++) mat[nfxs*ii+j] = 0;" << endl;
                    ofs << "for(int i=0; i<nfxs;i++) mat[nfxs*i+ii] = 0;" << endl;
                    ofs << "mat[ii*nfxs+ii] = 1.L-(1.L-2.L*x[ii])*r[ii];" << endl;
                }
                ofs  << "det = MatDetL(mat, nfxs);" << endl;
                
                ex intg = kv.second;
                cseParser cse;
                intg = cse.Parse(intg);
                ofs << "dCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                for(auto kv : cse.os()) {
                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                    Evalf(kv.second.subs(czRepl).subs(plRepl)).print(cppL);
                    ofs << ";" << endl;
                }
                
                ofs << "ytmp = ";
                Evalf(intg.subs(czRepl).subs(plRepl)).print(cppL);
                ofs << ";" << endl;
                ofs << "yy += det * ytmp;" << endl << endl;
                ofs << "}" << endl;
            }
        } else {
            auto tmp = expr.subs(FTX(wild(1),wild(2))==1).subs(cxRepl).subs(plRepl);
            if(SD::debug) {
                ofs << "//debug-int: " << tmp << endl;
            }
            ofs << "dCOMPLEX yy = ";
            Evalf(tmp).print(cppL);
            ofs << ";" << endl;
        }
        ofs << "y[0] = yy.real();" << endl;
        ofs << "y[1] = yy.imag();" << endl;
        ofs << "return 0;" << endl;
        ofs << "}" << endl;
/*----------------------------------------------*/
ofs << R"EOF(
#undef Pi
#undef Euler
#undef iEpsilon
#define Pi 3.1415926535897932384626433832795028841971693993751Q
#define Euler 0.57721566490153286060651209008240243104215933593992Q
#define iEpsilon 1.E-50Qi
)EOF" << endl;
/*----------------------------------------------*/
        // Quadruple
        auto cppQ = CppFormat(ofs, "Q");
        ofs << "extern \"C\" " << endl;
        ofs << "int SDQ_"<<rid<<"(const unsigned int xn, const qREAL x[], const int unsigned yn, qREAL y[], const qREAL pl[], const qREAL las[]) {" << endl;
        if(hasF) {
            ofs << "qREAL x0[xn];" << endl;
            ofs << "qCOMPLEX z[xn],r[xn];" << endl;
            ofs << "qREAL dff[xn+1];";
            ofs << "qCOMPLEX yy=0, ytmp, det;;" << endl;
            ofs << "int ii, nfxs="<<fxs.size()<<";" << endl;
            ofs << "qCOMPLEX mat[nfxs*nfxs];" << endl;
            for(auto &kv : ft_expr) {
                ofs << "{" << endl;
                lst xs0;
                for(int ii=0; ii<fxs.size(); ii++) {
                    if(!kv.first.has(fxs[ii])) xs0.append(ii);
                }
                ofs << "for(int i=0; i<xn; i++) z[i] = x0[i] = x[i];" << endl;
                for(auto x0i : xs0) ofs << "x0["<<x0i<<"]=0;" << endl;
                ofs << "X2ZQ_"<<ft_n<<"(x0,z,r,dff,pl,las);" << endl;
                ofs << "MatQ_"<<ft_n<<"(mat,x0,dff,pl,las);" << endl;
                for(auto x0i : xs0) {
                    ofs << "ii = " << x0i << ";" << endl;
                    ofs << "z[ii] = x[ii]-x[ii]*(1-x[ii])*r[ii];" << endl;
                    ofs << "for(int j=0; j<nfxs;j++) mat[nfxs*ii+j] = 0;" << endl;
                    ofs << "for(int i=0; i<nfxs;i++) mat[nfxs*i+ii] = 0;" << endl;
                    ofs << "mat[ii*nfxs+ii] = 1.Q-(1.Q-2.Q*x[ii])*r[ii];" << endl;
                }
                ofs  << "det = MatDetQ(mat, nfxs);" << endl;
                
                ex intg = kv.second;
                cseParser cse;
                intg = cse.Parse(intg);
                ofs << "qCOMPLEX "<<cse.oc<<"[" << cse.on()+1 << "];" << endl;
                for(auto kv : cse.os()) {
                    ofs <<cse.oc<< "["<<kv.first<<"] = ";
                    Evalf(kv.second.subs(czRepl).subs(plRepl)).print(cppQ);
                    ofs << ";" << endl;
                }
                
                ofs << "ytmp = ";
                Evalf(intg.subs(czRepl).subs(plRepl)).print(cppQ);
                ofs << ";" << endl;
                ofs << "yy += det * ytmp;" << endl << endl;
                ofs << "}" << endl;
            }
        } else {
            auto tmp = expr.subs(FTX(wild(1),wild(2))==1).subs(cxRepl).subs(plRepl);
            ofs << "qCOMPLEX yy = ";
            Evalf(tmp).print(cppQ);
            ofs << ";" << endl;
        }
        ofs << "y[0] = crealq(yy);" << endl;
        ofs << "y[1] = cimagq(yy);" << endl;
        ofs << "return 0;" << endl;
        ofs << "}" << endl;
        ofs << endl;

        ofs.close();
        
        ostringstream ofn, cmd;
        ofn << pid << "/" << rid << ".o";
        cmd << "g++ -fPIC " << CFLAGS << " -c -o " << ofn.str() << " " << cppfn.str();
        system(cmd.str().c_str());
        if(!debug) remove(cppfn.str().c_str());
        return lst{ rid, xs.size(), kvf.op(0), ft_n };
    }, "CI-I", Verbose, false);
    

//============================================================================================================

    if(use_cpp) {
        ostringstream sofn, garfn, cmd;
        if(key != NULL) {
            sofn << key << ".so";
            garfn << key << ".ci.gar";
            lst gar_res;
            for(auto &item : res) gar_res.append(item);
            archive ar;
            ar.archive_ex(gar_res, "res");
            ar.archive_ex(19790923, "c");
            ar.archive_ex(FT_N_NX, "ftnxn");
            ofstream out(garfn.str());
            out << ar;
            out.close();
        } else {
            sofn << pid << ".so";
            for(auto &item : res) ciResult.push_back(ex_to<lst>(item));
        }
        
        CompileMatDet();
        cmd << "g++ -rdynamic -fPIC -shared -lquadmath " << CFLAGS << " -o " << sofn.str() << " " << pid << "/*.o";
        system(cmd.str().c_str());
        cmd.clear();
        cmd.str("");
        cmd << "rm -rf " << pid;
        if(!debug) system(cmd.str().c_str());
    } else {
        // Using GiNaC internally
        ostringstream garfn;
        lst fn_lst;
        for(auto &kv : ftnvec) fn_lst.append(lst{kv.first, kv.second});
        
        //deleted from GiNaC 1.7.7
        //if(fn_lst.nops()<1) fn_lst = lst{lst{1, -1}};
        
        lst cifn_lst;
        for(auto &item : res_vec) cifn_lst.append(item);
        if(key != NULL) {
            garfn << key << ".ci.gar";
            archive ar;
            ar.archive_ex(fn_lst, "fn_lst");
            ar.archive_ex(cifn_lst, "cifn_lst");
            ar.archive_ex(FT_N_NX, "ftnxn");
            ar.archive_ex(19790923, "c");
            ofstream out(garfn.str());
            out << ar;
            out.close();
        } else {
            ciResult.clear();
            ciResult.shrink_to_fit();
            ciResult.push_back(fn_lst);
            ciResult.push_back(cifn_lst);
        }
    }
}

// need Parameter
void SD::Contours(const char *key, const char *pkey) {
    if(IsZero) return;
    
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    if(key != NULL) {
        ostringstream garfn;
        garfn << key << ".ci.gar";
        archive ar;
        ifstream in(garfn.str());
        in >> ar;
        in.close();
        auto c = ar.unarchive_ex(ParallelSymbols, "c");
        if(c!=19790923) throw runtime_error("*.ci.gar error!");
        FT_N_NX = ex_to<lst>(ar.unarchive_ex(ParallelSymbols, "ftnxn"));
        if(!use_cpp) {
            lst fn_lst = ex_to<lst>(ar.unarchive_ex(ParallelSymbols, "fn_lst"));
            ciResult.clear();
            ciResult.push_back(fn_lst);
        }
    }
    
    //change 2->1 from GiNaC 1.7.7
    if(FT_N_NX.nops()<1) return;
    if(Verbose > 0) cout << now() << " - Contours ..." << endl << flush;
    
    if(!use_cpp) {
        lst fn_lst = ciResult[0];
        exmap fn_map;
        for(auto item : fn_lst) fn_map[item.op(1)] = item.op(0);
        //change 1->0 from GiNaC 1.7.7
        for(int i=0; i<FT_N_NX.nops(); i++) FT_N_NX.let_op(i)= lst{FT_N_NX.op(i).op(0), fn_map[FT_N_NX.op(i).op(0)]};
        //now FT_N_NX item: lst { ft-id, ft}
    }
    
    vector<ex> nxn_vec;
    //change 1->0 from GiNaC 1.7.7
    for(int i=0; i<FT_N_NX.nops(); i++) nxn_vec.push_back(FT_N_NX.op(i));
    
    auto pid = getpid();
    ostringstream cmd;
    cmd << "mkdir -p " << pid;
    system(cmd.str().c_str());
    
    vector<ex> res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, nxn_vec, [&](auto & nxn, auto rid) {
        // return lst{ ft_n, lst{lambda-i, lambda-max} }
        // with I*[lambda-i]*lambda, lambda < lambda-max
        // note that lambda sequence only matches to x sequence in F-term
        int npara = -1;
        for(auto kv : Parameter) if(npara<kv.first) npara = kv.first;
        dREAL paras[npara+1];
        for(auto kv : Parameter) paras[kv.first] = CppFormat::ex2q(kv.second);
        
        int nvars = 0;
        ostringstream fname;
        void* module = nullptr;
        ex ft = 0; vector<ex> xs;
        if(use_cpp) {
            nvars = ex_to<numeric>(nxn.op(1)).to_int();
            ostringstream sofn;
            if(key != NULL) {
                sofn << key << ".so";
            } else {
                sofn << pid << ".so";
            }
            module = dlopen(sofn.str().c_str(), RTLD_NOW);
            if (module == nullptr) throw std::runtime_error("could not open compiled module!");
        } else {
            omp_set_num_threads(1);
            lst plReplace;
            for(auto kv : Parameter) plReplace.append(PL(kv.first)==kv.second);
            ft = nxn.op(1).subs(plReplace);
            xs = get_xy_from(ft);
            nvars = xs.size();
        }
        
        dREAL nlas[nvars];
        dREAL max_df = -1, max_f;
        for(int i=0; i<nvars; i++) {
            MinimizeBase::FunctionType dfp = NULL;
            if(use_cpp) {
                fname.clear();
                fname.str("");
                fname << "dirF_"<<nxn.op(0)<<"_"<<i;
                dfp = (MinimizeBase::FunctionType)dlsym(module, fname.str().c_str());
                assert(dfp!=NULL);
            } else {
                auto xi = get_xy_from(ft)[i];
                ex df = mma_diff(ft, xi);
                GWrapper::InitMinFunction(-abs(df), xs, 1);
                dfp = GWrapper::MinFunction;
            }
            dREAL maxdf = Minimizer->FindMinimum(nvars, dfp, paras);
            maxdf = -maxdf;
            nlas[i] = maxdf;
            if(max_df<maxdf) max_df = maxdf;
        }
        
//TODO: add other schema
//--------------------------------------------------
        for(int i=0; i<nvars; i++) {
            if(nlas[i] > 5E-3 * max_df) nlas[i] = 1/nlas[i];
            else nlas[i] = 1/max_df;
        }
        
        if(true) {
            dREAL nlas2 = 0;
            for(int i=0; i<nvars; i++) {
                nlas2 += nlas[i] * nlas[i];
            }
            nlas2 = sqrt(nlas2);
            for(int i=0; i<nvars; i++) {
                nlas[i] = nlas[i]/nlas2;
            }
        }
//--------------------------------------------------
        
        lst las;
        for(int i=0; i<nvars; i++) {
            las.append(CppFormat::q2ex(nlas[i]));
        }
        
        MinimizeBase::FunctionType fp;
        if(use_cpp) {
            fname.clear();
            fname.str("");
            fname << "imgF_"<<nxn.op(0);
            fp = (MinimizeBase::FunctionType)dlsym(module, fname.str().c_str());
            assert(fp!=NULL);
        } else {
            lst zRepl;
            for(int i=0; i<nvars; i++) {
                auto df = mma_diff(ft, xs[i]);
                ex zi = xs[i] - xs[i]*(1-xs[i])*df*I*CppFormat::q2ex(nlas[i])*x(nvars);
                zRepl.append(xs[i]==zi);
            }
            auto xs2 = xs;
            xs2.push_back(x(nvars));
            GWrapper::InitMinFunction(-ft.subs(zRepl), xs2, 2);
            fp = GWrapper::MinFunction;
        }
        
        dREAL laBegin = 0, laEnd = CTMax, min;
        dREAL UB[nvars+1];
        for(int i=0; i<nvars+1; i++) UB[i] = 1;
        
        min = laEnd;
        while(true) {
            UB[nvars] = min;
            dREAL res = Minimizer->FindMinimum(nvars+1, fp, paras, nlas, UB, NULL, NULL, true, CTTryPTS, CTSavePTS);
            if(res < -1E-5) laEnd = min;
            else laBegin = min;
            
            if(laEnd < 1E-10) {
                cout << RED << "too small lambda!" << RESET << endl;
                break;
            }
            
            if(laEnd-laBegin <= 5.E-2*laEnd) break;
            min = (laBegin + laEnd) / 2.0;
        }
        min = laBegin;
        
        if(use_dlclose) dlclose(module);
        
        las.append(CppFormat::q2ex(min));
        if(Verbose>5) {
            auto oDigits = Digits;
            Digits = 3;
            cout << "\r                                                    \r";
            cout << "     λ: " << las.evalf() << endl;
            Digits = oDigits;
        }
        return lst{ nxn.op(0), las };
    
    }, "las", Verbose, !debug);
    
    ostringstream garfn;
    if(key != NULL) {
        garfn << key;
        if(pkey != NULL) garfn << "-" << pkey;
        garfn << ".las.gar";
        lst gar_res;
        for(auto &item : res) gar_res.append(item);
        archive ar;
        ar.archive_ex(gar_res, "res");
        ar.archive_ex(19790923, "c");
        ofstream out(garfn.str());
        out << ar;
        out.close();
    } else {
        for(auto &item : res) LambdaMap[item.op(0)] = item.op(1);
    }
    
    cmd.clear();
    cmd.str("");
    cmd << "rm -rf " << pid;
    if(!debug) system(cmd.str().c_str());
}

// need Parameter
void SD::Integrates(const char *key, const char *pkey, int kid) {
    if(IsZero) return;
    
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    if(Verbose > 0) cout << now() << " - Integrates ..." << endl << flush;
    if(!use_cpp) omp_set_num_threads(1);
    
    lst lstRE;
    auto pid = getpid();
    ostringstream sofn, cmd;
    if(key == NULL) {
        sofn << pid << ".so";
        if(!use_cpp) {
            auto cifn_lst = ciResult[1];
            ciResult.clear();
            for(auto item : cifn_lst) ciResult.push_back(ex_to<lst>(item));
        }
    } else {
        sofn << key << ".so";
        ostringstream garfn;
        garfn << key << ".ci.gar";
        archive ar;
        ifstream in(garfn.str());
        in >> ar;
        in.close();
        auto c = ar.unarchive_ex(ParallelSymbols, "c");
        if(c!=19790923) throw runtime_error("*.ci.gar error!");
        if(use_cpp) {
            auto res = ar.unarchive_ex(ParallelSymbols, "res");
            ciResult.clear();
            for(auto item : ex_to<lst>(res)) ciResult.push_back(ex_to<lst>(item));
        } else {
            auto res = ex_to<lst>(ar.unarchive_ex(ParallelSymbols, "cifn_lst"));
            ciResult.clear();
            for(auto item : res) ciResult.push_back(ex_to<lst>(item));
        }
        
        garfn.clear();
        garfn.str("");
        garfn << key;
        if(pkey != NULL) garfn << "-" << pkey;
        garfn << ".las.gar";
        if(file_exists(garfn.str().c_str())) {
            archive la_ar;
            ifstream la_in(garfn.str());
            la_in >> la_ar;
            la_in.close();
            auto la_c = la_ar.unarchive_ex(ParallelSymbols, "c");
            auto la_res = la_ar.unarchive_ex(ParallelSymbols, "res");
            if(la_c!=19790923) throw runtime_error("*.ci.gar error!");
            for(auto item : ex_to<lst>(la_res)) {
                LambdaMap[item.op(0)] = item.op(1);
            }
        }
        
        if(kid>0) {
            garfn.clear();
            garfn.str("");
            garfn << key;
            if(pkey != NULL) garfn << "-" << pkey;
            garfn << ".res.gar";
            assert(file_exists(garfn.str().c_str()));
            
            archive res_ar;
            ifstream res_in(garfn.str());
            res_in >> res_ar;
            res_in.close();
            auto res_c = res_ar.unarchive_ex(ParallelSymbols, "c");
            auto relst = res_ar.unarchive_ex(ParallelSymbols, "relst");
            if(res_c!=19790923) throw runtime_error("*.res.gar error with kid!");
            lstRE = ex_to<lst>(relst);
        }
    }
    
    void* module = nullptr;
    if(use_cpp) {
        module = dlopen(sofn.str().c_str(), RTLD_NOW);
        if (module == nullptr) throw std::runtime_error("could not open compiled module!");
        if(!debug && key == NULL) remove(sofn.str().c_str());
    }
    
    int npara = 0;
    lst plRepl;
    for(auto kv : Parameter) {
        plRepl.append(PL(kv.first)==kv.second);
        if(kv.first>npara) npara = kv.first;
    }
    plRepl.sort();
    plRepl.unique();
    
    int total = ciResult.size(), current = 0;
    qREAL stot = sqrtq(total*1.Q);
    ResultError = 0;
    for(auto &item : ciResult) {
        current++;
        if(kid>0 && current != kid) continue;
        if(Verbose > 1) {
            cout << "\r  \\--Evaluating [" <<current<<"/"<<total<< "] ... " << flush;
        }
        
        unsigned int xsize = 0;
        ex co, exint, exft;
        vector<ex> xs, fxs;
        if(use_cpp) {
            xsize = ex_to<numeric>(item.op(1)).to_int();
            if(xsize<1) {
                Digits = 35;
                if(kid>0) {
                    lstRE.let_op(kid-1) = VE(item.op(0).subs(plRepl).evalf(),0) * item.op(2).subs(plRepl);
                    break;
                }
                ResultError +=  VE(item.op(0).subs(plRepl).evalf(),0) * item.op(2).subs(plRepl);
                lstRE.append(VE(item.op(0).subs(plRepl).evalf(),0) * item.op(2).subs(plRepl));
                continue;
            }
            co = item.op(2).subs(plRepl).subs(iEpsilon==0);
        } else {
            co = item.op(0).subs(plRepl);
            exint = item.op(1).subs(plRepl);
            exft = item.op(2).subs(plRepl);
            xs = get_xy_from(exint);
            fxs = get_xy_from(exft);
            xsize = xs.size();
            
            if(xsize<1) {
                Digits = 35;
                if(kid>0) {
                    lstRE.let_op(kid-1) = VE(exint.evalf(),0) * co;
                    break;
                }
                ResultError +=  VE(exint.evalf(),0) * co;
                lstRE.append(VE(exint.evalf(),0) * co);
                continue;
            }
        }
        
        if(co.is_zero()) continue;
        assert(!co.has(PL(wild())));
        qREAL cmax = -1;
        int reim = 0;
        if(ReIm==3) reim = 3;
        co = mma_collect(co, eps);
        for(int si=co.ldegree(eps); si<=co.degree(eps); si++) {
            auto tmp = co.coeff(eps, si);
            assert(!tmp.has(eps));
            tmp = mma_collect(tmp, ep);
            for(int i=tmp.ldegree(ep); i<=tmp.degree(ep); i++) {
                Digits = 50;
                auto ccRes = tmp.coeff(ep, i).evalf().expand();
                lst css;
                css.append(ccRes);
                if(is_a<add>(ccRes)) {
                    for(auto item : ccRes) css.append(item);
                }
                
                for(int ci=0; ci<css.nops(); ci++) {
                    auto nt = css.op(ci).subs(log(vs)==1).subs(vs==1).subs(nReplacements).evalf();
                    if(!is_a<numeric>(nt)) {
                        cerr << "nt: " << nt << endl;
                        assert(false);
                    }
                    assert(!nt.has(ep));
                    if(ReIm!=3 && reim!=3) {
                        if(ex_to<numeric>(nt).imag()==0) {
                            if(reim==2) reim = 3;
                            else reim = 1;
                        } else if(ex_to<numeric>(nt).real()==0) {
                            if(reim==1) reim = 3;
                            else reim = 2;
                        } else {
                            reim = 3;
                        }
                    }
                    nt = abs(nt).evalf(); // no PL here, then nReplacements
                    
                    qREAL qnt = CppFormat::ex2q(nt);
                    if(qnt > cmax) cmax = qnt;
                }
            }
        }
        if(cmax<=0) {
            cerr << "cmax<=0 with co = " << co <<endl;
            assert(false);
        }
        if(reim!=3 && ReIm!=3) {
            if(reim==1 && ReIm==2) reim=2;
            else if(reim==2 && ReIm==2) reim=1;
        }
        
        if(Verbose > 3) cout << "XDim=" << xsize << ", EpsAbs=" << (double)(EpsAbs/cmax/stot) << "/" << (double)cmax << endl;
        
        auto las = LambdaMap[item.op(3)];
        if(las.is_zero() && item.op(3)>0) {
            cerr << "lambda with the key(ft_n=" << item.op(3) << ") is NOT found!" << endl;
            assert(false);
        }
        
        IntegratorBase::SD_Type fp = nullptr, fpQ = nullptr;
        if(use_cpp) {
            int rid = ex_to<numeric>(item.op(0)).to_int();
            ostringstream fname;
            fname << "SDD_" << rid;
            fp = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
            assert(fp!=NULL);
            fname.clear();
            fname.str("");
            fname << "SDQ_" << rid;
            fpQ = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
            assert(fpQ!=NULL);
        } else {
            ex intx = 0;
            if(is_a<lst>(las)) {
                exset ftxset;
                exint.find(FTX(wild(1),wild(2)), ftxset);
                lst ftxlst;
                for(auto it : ftxset) ftxlst.append(it);
                auto expr = mma_collect(exint, FTX(wild(1),wild(2)));
                lst ftx_lst;
                for(auto item : fxs) ftx_lst.append(item);
                for(auto item : ftxlst) {
                    ex cftx = item.op(1);
                    ex cc = exint.coeff(item);
                    lst x0Repl;
                    for(auto xi : cftx) {
                        if(!ftx_lst.has(xi)) x0Repl.append(xi==0);
                    }
                    lst zRepl;
                    matrix mat(fxs.size(), fxs.size());
                    for(int i=0; i<fxs.size(); i++) {
                        ex df = mma_diff(exft, fxs[i]).subs(x0Repl);
                        ex zi = fxs[i] - fxs[i]*(1-fxs[i])*df*las.op(i)*I*z(0);
                        zRepl.append(fxs[i]==zi);
                        for(int j=0; j<fxs.size(); j++) mat(i, j) = mma_diff(zi, fxs[j]);
                    }
                    intx += mat.determinant() * cc.subs(zRepl);
                }
            } else {
                intx = exint.subs(FTX(wild(1),wild(2))==1);
            }
            GWrapper::InitIntFunction(intx, xs);
            fp = fpQ = GWrapper::IntFunction;
        }
        
        qREAL lambda[las.nops()];
        qREAL paras[npara+1];
        for(auto kv : Parameter) paras[kv.first] = CppFormat::ex2q(kv.second);
        
        Integrator->Verbose = Verbose;
        Integrator->ReIm = reim;
        if(is_a<lst>(las)) {
            qREAL lamax = CppFormat::ex2q(las.op(las.nops()-1));
            if(lamax > LambdaMax) lamax = LambdaMax;
            
            if(TryPTS<10000) TryPTS = 10000;
            Integrator->RunMAX = -100;
            Integrator->RunPTS = TryPTS/100;
            Integrator->EpsAbs = EpsAbs/cmax/2/stot;
            Integrator->EpsRel = 0;
            
            int smin = -1;
            int ctryR = 0, ctry = 0, ctryL = 0;
            ex cerr;
            qREAL log_lamax = log10q(lamax);
            qREAL log_lamin = log_lamax-1.Q;
            
            ostringstream las_fn;
            las_fn << key;
            if(pkey != NULL) las_fn << "-" << pkey;
            las_fn << "-" << current << ".las";
            if(use_las && file_exists(las_fn.str().c_str())) {
                std::ifstream las_ifs;
                las_ifs.open(las_fn.str(), ios::in);
                if (!las_ifs) throw runtime_error("failed to open *.las file!");
                for(int i=0; i<las.nops()-1; i++) {
                    dREAL la_tmp;
                    las_ifs >> la_tmp;
                    lambda[i] = la_tmp;
                }
                las_ifs.close();
                auto res = Integrator->Integrate(xsize, fp, fpQ, paras, lambda);
                auto res_tmp = res.subs(VE(wild(1), wild(2))==wild(2));
                auto err = real_part(res_tmp);
                if(err < imag_part(res_tmp)) err = imag_part(res_tmp);
                ErrMin::lastResErr = res;
                cerr = err;
            } else {
            // ---------------------------------------
            while(true) {
                smin = -1;
                ex emin = 0;
                ex lastResErr = 0;
                for(int s=0; s<=LambdaSplit; s++) {
                    if(Verbose>10 && s==0) {
                        if(ctryR>0 || ctry>0 || ctryL>0)
                            cout << "     ------------------------------" << endl;
                    }
                    auto log_cla = (log_lamin + s * (log_lamax-log_lamin) / LambdaSplit);
                    auto cla = powq(10.Q, log_cla);
                    if(cla < 1E-10) continue;
                    for(int i=0; i<las.nops()-1; i++) {
                        lambda[i] = CppFormat::ex2q(las.op(i)) * cla;
                    }
 
                    if(!use_cpp) GWrapper::Lambda = CppFormat::q2ex(cla);
                    auto res = Integrator->Integrate(xsize, fp, fpQ, paras, lambda);
                    if(Verbose>10) {
                        auto oDigits = Digits;
                        Digits = 3;
                        cout << "\r                                                    \r";
                        if(res.has(NaN)) cout << "     λ=" << (double)cla << ": " << NaN << endl;
                        else cout << "     λ=" << (double)cla << ": " << HepLib::VEResult(VESimplify(res)) << endl;
                        Digits = oDigits;
                    }
                    
                    if(res.has(NaN) && s==0) continue;
                    else if(res.has(NaN)) break;
                    if(lastResErr.is_zero()) lastResErr = res;
                    auto diff = VESimplify(lastResErr - res);
                    diff = diff.subs(VE(0,0)==0);
                    exset ves;
                    diff.find(VE(wild(0), wild(1)), ves);
                    bool err_break = false;
                    for(auto ve : ves) {
                        if(abs(ve.op(0)) > 3*ve.op(1)) {
                            err_break = true;
                            break;
                        }
                    }
                    if(err_break) break;
                    lastResErr = res;
                    
                    auto res_tmp = res.subs(VE(wild(1), wild(2))==wild(2));
                    auto err = real_part(res_tmp);
                    if(err < imag_part(res_tmp)) err = imag_part(res_tmp);
                    if(smin<0 || err < emin) {
                        ErrMin::lastResErr = res;
                        cerr = err;
                        smin = s;
                        emin = err;
                        if(emin < CppFormat::q2ex(EpsAbs/cmax/stot)) {
                            smin = -2;
                            if(kid>0) {
                                lstRE.let_op(kid-1) = co * res;
                                break;
                            }
                            ResultError += co * res;
                            lstRE.append(co * res);
                            if(Verbose>5) {
                                cout << WHITE;
                                cout << "     λ=" << (double)cla << ": " << HepLib::VEResult(VESimplify(res)) << endl;
                                cout << RESET;
                            }
                            break;
                        }
                    } else if(err > 1.E3 * emin) {
                        break;
                    }
                }
                if(smin == -2) break;
                if(smin == -1) {
                    std::cerr << "smin = -1, optimized lambda NOT found!" << endl;
                    assert(false);
                }
                
                if(smin <= 0) {
                    if(ctryL >= CTryLeft || ctryR>0) break;
                    log_lamax = log_lamin;
                    log_lamin -= 1.5Q;
                    ctryL++;
                } else if(smin >= LambdaSplit) {
                    if(ctryR >= CTryRight || ctryL>0) break;
                    log_lamin = log_lamax;
                    log_lamax += log10q(CTryRightRatio);
                    ctryR++;
                } else {
                    if(ctry >= CTry) break;
                    auto la1 = log_lamin + (smin-1) * (log_lamax-log_lamin) / LambdaSplit;
                    auto la2 = log_lamin + (smin+1) * (log_lamax-log_lamin) / LambdaSplit;
                    log_lamin = la1;
                    log_lamax = la2;
                    ctry++;
                }
            }
            
            if(smin == -2) {
                if(kid>0) break;
                continue;
            }
            
            auto log_cla = (log_lamin + smin * (log_lamax-log_lamin) / LambdaSplit);
            auto cla = powq(10.Q, log_cla);
            if(Verbose > 7) cout << WHITE << "     Final λ = " << (double)cla << " / " << las.op(las.nops()-1) << RESET << endl;
            for(int i=0; i<las.nops()-1; i++) {
                lambda[i] = CppFormat::ex2q(las.op(i)) * cla;
            }
            // ---------------------------------------
            }
            
            // ---------------------------------------
            // try HookeJeeves
            // ---------------------------------------
            if( use_ErrMin && (cerr > CppFormat::q2ex(1E5 * EpsAbs/cmax/stot)) ) {
                auto miner = new HookeJeeves();
                ErrMin::miner = miner;
                ErrMin::xsize = xsize;
                ErrMin::fp = fp;
                ErrMin::fpQ = fpQ;
                ErrMin::paras = paras;
                ErrMin::Integrator = Integrator;
                ErrMin::Verbose = Verbose;
                dREAL oo[las.nops()-1], ip[las.nops()-1];
                for(int i=0; i<las.nops()-1; i++) ip[i] = oo[i] = lambda[i];
                ErrMin::lambda = oo;
                ErrMin::err_max = 1E30;
                auto err_min = ErrMin::err_min;
                ErrMin::err_min = err_min < 0 ? -err_min * ex_to<numeric>(cerr).to_double() : err_min/cmax;
                ErrMin::RunRND = 0;
                miner->Minimize(las.nops()-1, ErrMin::IntError, ip);
                delete miner;
                ErrMin::err_min = err_min;
                for(int i=0; i<las.nops()-1; i++) {
                    lambda[i] = ErrMin::lambda[i];
                }
            }
            // ---------------------------------------
            
            if(save_las) {
                std::ofstream las_ofs;
                las_ofs.open(las_fn.str(), ios::out);
                if (las_ofs) {
                    for(int i=0; i<las.nops()-1; i++) {
                        dREAL la_tmp = lambda[i];
                        las_ofs << la_tmp << " ";
                    }
                    las_ofs << endl;
                    las_ofs.close();
                }
            }
        }

        Integrator->RunMAX = RunMAX;
        Integrator->RunPTS = RunPTS;
        Integrator->EpsAbs = EpsAbs/cmax/stot;
        Integrator->EpsRel = 0;
        auto res = Integrator->Integrate(xsize, fp, fpQ, paras, lambda);
        if(Verbose>5) { 
            cout << WHITE;
            cout << "     Res = "<< HepLib::VEResult(VESimplify(res)) << endl;
            cout << RESET;
        }
        if(res.has(NaN)) {
            ResultError = NaN;
            if(kid>0) {
                lstRE.let_op(kid-1) = NaN;
            } else {
                lstRE.append(NaN);
            }
            break;
        } else {
            if(kid>0) {
                lstRE.let_op(kid-1) = co * res;
                break;
            } else {
                ResultError += co * res;
                lstRE.append(co * res);
            }
        }
    }
    
    if(use_dlclose) dlclose(module);
    if(total>0 && Verbose > 1) cout << "@" << now(false) << endl;
    
    if(kid>0) {
        ResultError = 0;
        for(auto item : lstRE) ResultError += item;
    }
    ResultError = VESimplify(ResultError,epN,epsN);
    
    if(key != NULL) {
        ostringstream garfn;
        garfn << key;
        if(pkey != NULL) garfn << "-" << pkey;
        garfn << ".res.gar";
        archive ar;
        ar.archive_ex(ResultError, "res");
        ar.archive_ex(lstRE, "relst");
        ar.archive_ex(19790923, "c");
        ofstream out(garfn.str());
        out << ar;
        out.close();
    }
}

void SD::Initialize(XIntegrand xint) {
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    IsZero = false;
    Replacements2(xint.nReplacements);
    
    nReplacements = xint.nReplacements;
    Deltas = xint.Deltas;

    FunExp.clear();
    FunExp.push_back(make_pair(xint.Functions, xint.Exponents));
    
    // Do Other Simplifications
    Normalizes();
    XReOrders();
    Normalizes();
}

void SD::Evaluate(FeynmanParameter fp, const char* key) {
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    cout << endl << "Starting @ " << now() << endl;
    if(SecDec==NULL) SecDec = new SecDecG();
    if(Integrator==NULL) Integrator = new HCubature();
    if(Minimizer==NULL) Minimizer = new MinUit();
    
    Initialize(fp);
    if(FunExp.size()<1) return;
    Scalelesses();
        
    RemoveDeltas();
    SDPrepares();
    EpsEpExpands();
    CIPrepares(key);
    Contours(key);
    Integrates(key);
    delete SecDec;
    delete Integrator;
    delete Minimizer;
    cout << "Finished @ " << now() << endl << endl;
}

void SD::Evaluate(XIntegrand xint, const char *key) {
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    cout << endl << "Starting @ " << now() << endl;
    if(SecDec==NULL) SecDec = new SecDecG();
    if(Integrator==NULL) Integrator = new HCubature();
    
    //if(Minimizer==NULL) Minimizer = new HookeJeeves();
    if(Minimizer==NULL) Minimizer = new MinUit();
    
    Initialize(xint);
    if(FunExp.size()<1) return;
    
    RemoveDeltas();
    SDPrepares();
    EpsEpExpands();
    CIPrepares(key);
    Contours(key);
    Integrates(key);
    delete SecDec;
    delete Integrator;
    delete Minimizer;
    cout << "Finished @ " << now() << endl << endl;
}

double SD::FindMinimum(ex expr, bool compare0) {
    static long long fid = 0;
    fid++;
    ostringstream cppfn, sofn, cmd;
    auto pid = getpid();
    cppfn << "/tmp/" << pid << "-" << fid << "-min.cpp";
    sofn << "/tmp/" << pid << "-" << fid << "-min.so";
    std::ofstream ofs;
    ofs.open(cppfn.str(), ios::out);
    if (!ofs) throw runtime_error("failed to open *.cpp file!");
    
    auto xs = get_xy_from(expr);
    dREAL UB[xs.size()+ParameterUB.size()], LB[xs.size()+ParameterUB.size()];
    lst cxRepl;
    int count = 0;
    for (auto xi : xs) {
        ostringstream xs;
        xs << "x[" << count << "]";
        cxRepl.append(xi == symbol(xs.str()));
        UB[count] = 1;
        LB[count] = 0;
        count++;
    }
    for(auto kv : ParameterUB) {
        int key = kv.first;
        ostringstream xs;
        xs << "x[" << count << "]";
        cxRepl.append(PL(key) == symbol(xs.str()));
        UB[count] = ParameterUB[key].to_double();
        LB[count] = ParameterLB[key].to_double();
        count++;
    }
    
/*----------------------------------------------*/
ofs << R"EOF(
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <iostream>
using namespace std;

typedef long double dREAL;
typedef complex<long double> dCOMPLEX;

#define Pi 3.1415926535897932384626433832795028841971693993751L
#define Euler 0.57721566490153286060651209008240243104215933593992L

//#define expt(a,b) pow(a,b)
//#define recip(a) pow(a,-1)
dREAL expt(dREAL a, dREAL b) { return pow(a,b); }
dCOMPLEX expt(dCOMPLEX a, dREAL b) { return pow(a,b); }
dREAL recip(dREAL a) { return 1.L/a; }
dCOMPLEX recip(dCOMPLEX a) { return 1.L/a; }

)EOF" << endl;
/*----------------------------------------------*/

    auto cppL = CppFormat(ofs, "L");
    ofs << "extern \"C\" " << endl;
    ofs << "dREAL minFunc(int xn, dREAL* x, dREAL *pl, dREAL *las) {" << endl;
    auto tmp = expr.subs(cxRepl);
    assert(!tmp.has(PL(wild())));
    ofs << "dREAL yy = ";
    tmp.print(cppL);
    ofs << ";" << endl;
    ofs << "return yy;" << endl;
    ofs << "}" << endl;
    
    cmd.clear();
    cmd.str("");
    cmd << "g++ -fPIC -shared " << CFLAGS << " -o " << sofn.str() << " " << cppfn.str();
    system(cmd.str().c_str());
    
    void* module = nullptr;
    module = dlopen(sofn.str().c_str(), RTLD_NOW);
    if(module == nullptr) throw std::runtime_error("could not open compiled module!");
    
    auto fp = (MinimizeBase::FunctionType)dlsym(module, "minFunc");
    assert(fp!=NULL);
    
    double min = Minimizer->FindMinimum(count, fp, NULL, NULL, UB, LB, NULL, compare0, 5, 5);
    
    if(use_dlclose) dlclose(module);
    cmd.clear();
    cmd.str("");
    cmd << "rm " << cppfn.str() << " " << sofn.str();
    system(cmd.str().c_str());
    return min;
}

int SD::epRank(ex expr_in) {
    if(!expr_in.has(ep)) return 0;
    int p = -5;
    auto expr = mma_collect(expr_in, ep);
    while(true) {
        auto tmp = normal(series_to_poly(expr.series(ep, p)));
        if(!tmp.is_zero()) {
            tmp = mma_collect(tmp, ep);
            return tmp.ldegree(ep);
        } else p++;
    }
    throw runtime_error("epRank error!");
}

int SD::epsRank(ex expr_in) {
    if(!expr_in.has(eps)) return 0;
    int p = -5;
    auto expr = mma_collect(expr_in, eps);
    while(true) {
        auto tmp = normal(series_to_poly(expr.series(eps, p)));
        if(!tmp.is_zero()) {
            tmp = mma_collect(tmp, eps);
            return tmp.ldegree(eps);
        } else p++;
    }
    throw runtime_error("epsRank error!");
}

ex SD::VEResult() {
    return HepLib::VEResult(ResultError);
}

void SD::VEPrint(bool endlQ) {
    ex expr = HepLib::VEResult(ResultError);
    for(int i=expr.ldegree(eps); i<=expr.degree(eps); i++) {
        ex exp1 = expr.coeff(eps, i);
        for(int j=expr.ldegree(ep); j<=expr.degree(ep); j++) {
            cout << WHITE <<"(" << RESET;
            cout << exp1.coeff(ep, j);
            cout << WHITE << ")" << RESET;
            if(j!=0 || i!=0) cout << "*" << WHITE << pow(ep,j)*pow(eps,i) << RESET;
            if(j<expr.degree(ep)) cout << " + ";
        }
    }
    if(endlQ) cout << endl;
}

ex SD::PrefactorFIESTA(int nLoop) {
    return  pow(I*pow(Pi,2-ep)*exp(-ep*Euler), -ex(nLoop));
}

}
