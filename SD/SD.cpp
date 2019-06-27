#include "SD.h"
#include <math.h>

namespace HepLib {

symbol const SD::ep("ep");
symbol const SD::eps("eps");
symbol const SD::iEpsilon("iEpsilon");
realsymbol const SD::NaN("NaN");
bool SD::use_dlclose = true;

vector<exmap> SecDecBase::x2y(const lst &xpols) {
    ex xpol = 1;
    for(int i=0; i<xpols.nops(); i++) xpol *= xpols.op(i);
    return x2y(xpol);
}

/*********************************************************/
/*                    CheckFAtX1                         */
/*********************************************************/
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
        for(int pi=0; pi<std::pow(2, ys_tmp.size()); pi++) {
            lst yrepl;
            int ci = pi;
            for(int i=0; i<ys_tmp.size(); i++) {
                yrepl.append(ys_tmp[i] == ((ci % 2) ? 1 : 0));
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
    
    cout << "F: " << po_ex.first.op(1) << endl;
    cout << "F1 Failed with ALL possible bisections!" << endl;
    assert(false);
    return vector<pair<lst, lst>>();
}

/*********************************************************/
/*                    SDPrepare                         */
/*********************************************************/
// 1st element in input polist is the constant term
// 2nd element in input polist is the F-term
vector<pair<exmap, ex>> SD::SDPrepare(pair<lst, lst> po_ex) {
    lst const polist = po_ex.first;
    lst const exlist = po_ex.second;
    lst sdList;
    for(int i=0; i<polist.nops(); i++) {
        auto tmp = polist.op(i);
        auto ntmp = exlist.op(i);
        if(!tmp.subs(lst{x(wild())==0, y(wild())==0}).normal().is_zero()) continue;
        if( (!tmp.has(x(wild())) && !tmp.has(y(wild()))) || (is_a<numeric>(ntmp) && ntmp.evalf()>0) ) continue;
        sdList.append(tmp);
    }
    
    vector<exmap> vmap = SecDec->x2y(sdList);

    vector<pair<exmap, ex>> sd;
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

        bool need_contour_deformation = ft.has(PL(wild()));
        if(ft.has(y(wild())) && !need_contour_deformation) {
            auto tmp = ft.subs(nReplacements).expand();
            if(is_a<add>(tmp)) {
                need_contour_deformation = false;
                auto first = tmp.op(0).subs(y(wild())==1);
                for(auto item : tmp) {
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
                cout << "tmp = " << tmp << endl;
                cout << "tmp is NOT a numeric." << endl;
                assert(false);
            }
            if( tmp.evalf() < 0 ) {
                ct = exp(-I * Pi * exlist.op(1));
                fsgin = -1;
            }
            ft = 1;
        } else {
            auto tmp = ft.subs(nReplacements);
            auto tn = tmp.subs(y(wild())==numeric(0.5));
            for(auto kv : ParameterUB) {
                int k = kv.first;
                tn = tn.subs(PL(k)==(ParameterUB[k]+ParameterLB[k])/ex(2));
            }
            if(tn.has(PL(wild()))) {
                cout << "PL still exists in " << tn << endl;
                assert(false);
            }
            if(!is_a<numeric>(tn.evalf())) {
                cout << "tn is not numeric: " << tn.eval() << endl;
                assert(false);
            }
            
            if(tn.evalf() < 0) tmp = ex(0)-tmp;
            double tmin = FindMinimum(tmp);
            if(tmin > 0) {
                if(tn < 0) {
                    ct = exp(-I * Pi * exlist.op(1));
                    fsgin = -1;
                }
                ft = 1;
                need_contour_deformation = false;
            }
        }

        ex polexp = FT(ft)*CT(ct);
        exmap ymol;
        
        // need collect_common_factors
        auto det = collect_common_factors(vi[x(-1)]);
        assert(!is_exactly_a<add>(det));
        auto ys = get_xy_from(det);
        for(int i=0; i<ys.size(); i++) {
            ymol[ys[i]] = ymol[ys[i]] + det.degree(ys[i]);
        }

        for(int i=0; i<ypolist.nops(); i++) {
            // need collect_common_factors
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
            
            if(tmp.has(y(wild()))) tmp = collect_common_factors(tmp.expand());

            if(is_exactly_a<mul>(tmp)) {
                ex rem = 1;
                ex ct = 1;
                for (auto item : tmp) {
                    if( item.match(y(wild())) || item.match(pow(y(wild()), wild(1))) ) {
                        auto yi = get_xy_from(item)[0];
                        ymol[yi] = ymol[yi] + (item.nops()<2 ? 1 : item.op(1)) * exlist.op(i);
                    } else if(!item.has(y(wild())) && !item.has(x(wild()))) {
                        if(is_a<numeric>(exlist.op(i)) && ex_to<numeric>(exlist.op(i)).is_integer()) {
                            ct *= item;
                        } else if(!item.has(PL(wild()))) {
                            auto tr = item.subs(nReplacements);
                            if(!is_a<numeric>(tr.evalf())) {
                                cout << "not numeric - item: " << tr << " ; " << item << endl;
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
                        if(nchk && item.subs(y(wild())==0).subs(iEpsilon==0).is_zero()) {
                            cout << "zero - item: " << item << endl;
                            cout << "exlist.op(i) = " << exlist.op(i) << endl;
                            assert(false);
                        }
                        rem *= item;
                    }
                }
                polexp *= pow(rem-(i==1 && ft!=1 ? iEpsilon : ex(0)), exlist.op(i));
                polexp = polexp.subs(CT(wild()) == CT(wild()*pow(ct, exlist.op(i))));
            } else {
                auto item = tmp;
                if( item.match(y(wild())) || item.match(pow(y(wild()), wild(1))) ) {
                    auto yi = get_xy_from(item)[0];
                    ymol[yi] = ymol[yi] + (item.nops()<2 ? 1 : item.op(1)) * exlist.op(i);
                } else if(!tmp.has(x(wild())) && !tmp.has(y(wild()))) {
                    polexp = polexp.subs(CT(wild()) == CT(wild()*pow(tmp, exlist.op(i))));
                } else {
                    if(nchk && tmp.subs(y(wild())==0).is_zero()) {
                        cout << "zero - tmp: " << tmp << endl;
                        assert(false);
                    }
                    polexp *= pow(tmp-(i==1 && ft!=1 ? iEpsilon : ex(0)), exlist.op(i));
                }
            }
        }

        exmap xmol;
        for(auto & kv : ymol) {
            auto k = kv.first.subs({y(wild())==x(wild())});
            auto v = kv.second.subs({y(wild())==x(wild())});
            if(is_a<numeric>(v)) {
                auto nv = ex_to<numeric>(v);
                if(nv<=-1) {
                    cout << endl << "nv: " << nv << ", NOT larger than -1." << endl;
                    assert(false);
                }
            }
            xmol[k] = v;
        }
        polexp = polexp.subs({y(wild())==x(wild())});
        sd.push_back(make_pair(xmol, polexp));
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
                        } else if(!tmp.has(PL(wild()))) {
                            auto tr = tmp.subs(nReplacements);
                            if(!is_a<numeric>(tr)) {
                                cout << "tmp: " << tmp << endl;
                                cout << "tmp is NOT numeric with nReplacements." << endl;
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

void inline Replacements2(exmap &repl) {
    auto tmp = repl;
    for(auto &kv : repl) {
        kv.second = kv.second.subs(tmp, subs_options::algebraic);
    }
}

bool inline xPositive(ex expr) {
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
    if(fp.Propagators.nops() != fp.Exponents.nops()) {
        cout << "the length of Propagators and Exponents are NOT equal." << endl;
        assert(false);
    }
    
    if(fp.Prefactor.is_zero()) {
        IsZero = true;
        return;
    }
    
    IsZero = false;
    
    for(auto kv: fp.lReplacements) assert(!(lst{kv.first, kv.second}).has(iEpsilon));
    for(auto kv: fp.tReplacements) assert(!(lst{kv.first, kv.second}).has(iEpsilon));
    
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
                assert(is_a<numeric>(p.coeff(m,2)));
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
        ex usgn = u_nd.op(1).subs(xtNeg).subs(x(wild())==1/2).subs(nsubs).subs(plsubs);
        usgn = normal(usgn)>0 ? 1 : -1;
        
        if(!xPositive(normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacements).subs(plsubs))) {
            cout << "NOT positive - un: " << normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacements).subs(plsubs) << endl;
            assert(false);
        }
        if(!xPositive(normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacements).subs(plsubs))) {
            cout << "NOT positive - ud: " << normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacements).subs(plsubs) << endl;
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
        ex usgn = u_nd.op(1).subs(xtNeg).subs(x(wild())==1/2).subs(nsubs).subs(plsubs);
        usgn = normal(usgn)>0 ? 1 : -1;
        
        if(!xPositive(normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacements).subs(plsubs))) {
            cout << "NOT positive - un: " << normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacements).subs(plsubs) << endl;
            assert(false);
        }
        if(!xPositive(normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacements).subs(plsubs))) {
            cout << "NOT positive - ud: " << normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacements).subs(plsubs) << endl;
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
    
    ex usgn = u_nd.op(1).subs(xtNeg).subs(x(wild())==1/2).subs(nsubs).subs(plsubs);
    usgn = normal(usgn)>0 ? 1 : -1;
    ex fsgn = rem_nd.op(1).subs(xtNeg).subs(x(wild())==1/2).subs(nsubs).subs(plsubs);
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
                auto tmp = diff_wrt(kv, x(i));
                for(auto nkv : tmp) nret.push_back(nkv);
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

/*********************************************************/
/*               's Funtions in SD                       */
/*********************************************************/
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

void SD::Scalelesses() {
    if(IsZero) return;
    if(Deltas.size()<1) return;
    
    vector<pair<lst,lst>> ret;
    for(auto funexp : FunExp) {
        symbol s;
        auto fun = funexp.first;
        auto exp = funexp.second;
        bool is0;
        for(auto delta : Deltas) {
            is0 = false;
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
                    if(is_a<numeric>(exp.op(j)) && ex_to<numeric>(exp.op(j)).is_nonneg_integer() ) continue;
                    auto tmp = fun.op(j).subs(sRepl);
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
        if(!is0) ret.push_back(make_pair(fun, exp));
    }
    FunExp = ret;
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
                        cout << "xj: " << xj << endl;
                        cout << "fun: " << fun << endl;
                        cout << "fun is NOT polynormial of xj." << endl;
                        assert(false);
                    }
                    auto expn = fun.degree(xj);
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
    
    if(Verbose > 0) cout << now() << " - SDPrepares ..." << endl << flush;
    symbol s;
    Integrands.clear();
    
    if(CheckF1) {
        vector<ex> funexps =
        GiNaC_Parallel(ParallelProcess, ParallelSymbols, FunExp, [&](auto &kv, auto rid) {
            lst para_res_lst;
            auto kvs = AutoF1(kv);
            for(auto item : kvs) {
                para_res_lst.append(lst{item.first, item.second});
            }
            return para_res_lst;
        }, "f1", 0);

        FunExp.clear();
        for(auto &item : funexps) {
            for(auto &it : ex_to<lst>(item)) {
                FunExp.push_back(make_pair(ex_to<lst>(it.op(0)), ex_to<lst>(it.op(1))));
            }
        }
    }
    
    vector<ex> res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, FunExp, [&](auto &kv, auto rid) {

        // return single element in which ep/eps can be expanded safely.
        lst para_res_lst;
        auto xmol_exps = SDPrepare(kv);
        int run_count = 0;
                
        while(xmol_exps.size()>0) {
            if((++run_count) > 100) throw runtime_error("run count > 100 limit!");
            
            vector<pair<exmap, ex>> todo_list;
            for(auto const &xmol_exp : xmol_exps) {
                auto xmol = xmol_exp.first;
                auto expr = xmol_exp.second;
                bool pole_reached = true;
                
                exset fts;
                expr.find(FT(wild()), fts);
                bool noFT = (fts.size()==1) && ( (*(fts.begin())) == FT(1) );
                
                ex pole_requested = -1;
                if(noFT || PoleRequested > -1) pole_requested = PoleRequested;
                
                for(auto xen : xmol) {
                    auto expn = xen.second.subs(lst{eps==0,ep==0}).normal();
                    assert(is_a<numeric>(expn));
                    
                    if(ex_to<numeric>(expn) < pole_requested) {
                        pole_reached = false;
                        auto xx = xen.first;
                        expr = expr / (xen.second+1);
                        
                        exmap xmol2;
                        for(auto kv : xmol) {
                            if(kv.first != xx) xmol2[kv.first] = kv.second;
                        }
                        auto exp2 = expr.subs(xx==1);
                        todo_list.push_back(make_pair(xmol2, exp2));
                        
                        xmol[xx] = xen.second+1;
                        expr = ex(0)-expr.subs(xx==s).diff(s).subs(s==xx);
                        todo_list.push_back(make_pair(xmol, expr));
                        break;
                    }
                }
                
                if(pole_reached) {
                    lst exprs = {expr};
                    symbol dx;
                    for(auto xen : xmol) {
                        auto expn = xen.second.subs(lst{eps==0,ep==0}).normal();
                        assert(is_a<numeric>(expn));
                        
                        if(!noFT && ex_to<numeric>(expn)<-1) {
                            cout << "expn: " << expn << endl;
                            cout << "expn<-1 in the case with FT." << endl;
                            assert(false);
                        }

                        lst exprs2;
                        for(auto &it : exprs) {
                            ex rem = pow(xen.first, xen.second) * it;
                            if(ex_to<numeric>(expn)<=-1) {
                                ex dit = it;
                                ex dit0 = dit.subs(xen.first==0);
                                ex ifact = 1;
                                rem -= pow(xen.first, xen.second) * dit0 / ifact;
                                exprs2.append(dit0/(xen.second+1)/ifact);
                                for(int i=1; i+expn<0; i++) {
                                    dit = dit.subs(xen.first==dx).diff(dx).subs(dx==xen.first);
                                    dit0 = dit.subs(xen.first==0);
                                    ifact *= i;
                                    rem -= pow(xen.first, xen.second+i) * dit0 / ifact;
                                    exprs2.append(dit0/(xen.second+i+1)/ifact);
                                }
                            }
                            exprs2.append(rem);
                        }
                        exprs = exprs2;
                    }
                    
                    for(auto const &it : exprs) {
                        if(!it.is_zero()) para_res_lst.append(it);
                    }
                }
            }
            xmol_exps = todo_list;
        }
        
        for(int i=0; i< para_res_lst.nops(); i++) {
            exset xset;
            para_res_lst.op(i).find(x(wild()), xset);
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
            
            para_res_lst.let_op(i) = para_res_lst.op(i).subs(x2y).subs(y(wild())==x(wild()));
        }
        
        if(para_res_lst.nops()<1) para_res_lst.append(0);
        return para_res_lst;
    }, "sd", Verbose);
    
    for(auto &item : res) {
        for(auto &it : ex_to<lst>(item)) Integrands.push_back(it);
    }
}

void SD::EpsEpExpands() {
    if(IsZero) return;
    if(Integrands.size()<1) {
        IsZero = true;
        return;
    }
    
    int exN = 0;
    if(Verbose > 0) cout << now() << " - EpsEpExpands ..." << endl << flush;
    
    vector<ex> res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, Integrands, [&](auto &item, auto rid) {
        // return two elements,
        // 1st: x-independent coefficient, expanded in ep/eps
        // 2nd: x-integrand
        exset cts;
        item.find(CT(wild()), cts);
        if(cts.size() != 1) {
            cout << "CT size is NOT 1!" << endl;
            assert(false);
        }
        ex ct = (*(cts.begin())).subs(CT(wild())==wild());
        auto tmp = item.subs(CT(wild())==1);
        tmp = collect_common_factors(tmp); //TODO: need collect_common_factors ?
        lst para_res_lst;
        
        if(!tmp.has(eps) && !ct.has(eps)) {
            int ctN = epRank(ct);
            tmp = series_to_poly(tmp.series(ep, 1+exN+epN-ctN));
            for(int di=tmp.ldegree(ep); (di<=tmp.degree(ep) && di<=epN-ctN); di++) {
                auto intg = tmp.coeff(ep, di);
                auto pref = series_to_poly(ct.series(ep, 1+exN+epN-di));
                intg = collect_common_factors(intg); //TODO: collect_common_factors ?
                para_res_lst.append(lst{pref * pow(ep, di), intg});
            }
        } else {
            auto sct = ct;
            int sctN = epsRank(sct);
            ex stmp = series_to_poly(tmp.series(eps, 1+exN+epsN-sctN));
            stmp = collect_common_factors(stmp); //TODO: collect_common_factors ?
            for(int sdi=stmp.ldegree(eps); (sdi<=stmp.degree(eps) && sdi<=epsN-sctN); sdi++) {
                tmp = stmp.coeff(eps, sdi);
                tmp = collect_common_factors(tmp); //TODO: collect_common_factors ?
                assert(!tmp.has(eps));
                ct = series_to_poly(sct.series(eps, 1+exN+epsN-sdi));
                int ctN = epRank(ct);
                tmp = series_to_poly(tmp.series(ep, 1+exN+epN-ctN));
                for(int di=tmp.ldegree(ep); (di<=tmp.degree(ep) && di<=epN-ctN); di++) {
                    auto intg = tmp.coeff(ep, di);
                    assert(!intg.has(ep));
                    auto pref = series_to_poly(ct.series(ep, 1+exN+epN-di));
                    intg = collect_common_factors(intg); //TODO: collect_common_factors ?
                    para_res_lst.append(lst{pref * pow(eps, sdi) * pow(ep, di), intg});
                }
            }
        }
        
        if(para_res_lst.nops()<1) para_res_lst.append(lst{0,0});
        return para_res_lst;

    }, "ep", Verbose);
    
    epIntegrands.clear();
    for(auto &item : res) {
        for(auto &kv : ex_to<lst>(item)) {
            epIntegrands.push_back(make_pair(kv.op(0), kv.op(1)));
        }
    }
}

void SD::CompileMatDet() {
    auto pid = getpid();
    std::ofstream ofs;
    ostringstream cppfn, cmd;
    cppfn << pid << "/MatDet.cpp";
    ofs.open(cppfn.str(), ios::out);
    if (!ofs) throw runtime_error("failed to open *.cpp file!");
/****************************************************************/
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
dCOMPLEX MatDetD(dCOMPLEX mat[], int n) {
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
/****************************************************************/
    ofs.close();
    cmd.clear();
    cmd.str("");
    cmd << "g++ -fPIC " << CFLAGS << " -c -o " << pid << "/MatDet.o " << pid << "/MatDet.cpp";
    system(cmd.str().c_str());
}

void SD::IntPrepares(const char* key) {
    if(IsZero) return;
    if(epIntegrands.size()<1) {
        IsZero = true;
        return;
    }
    
    if(Verbose > 0) cout << now() << " - IntPrepares ..." << endl << flush;
    auto pid = getpid();
    
    vector<ex> res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, epIntegrands, [&](auto &kv, auto rid) {
        // return lst{ no-x-result, xn, x-indepent prefactor, F-term }
        // or     lst{ id(SD_D|Q[id] in .so), xn, x-indepent prefactor, F-term }
        auto expr = kv.second;
        auto xs = get_xy_from(expr);
        if(xs.size()<1) {
            return lst{
                expr.subs(FT(wild())==1).subs(iEpsilon==I*power(10,-30)),
                xs.size(), kv.first, 1
            };
        }

        lst cxRepl, czRepl, plRepl;
        int count = 0;
        for (auto xi : xs) {
            ostringstream xs, zs;
            xs << "x[" << count << "]";
            zs << "z[" << count << "]";
            cxRepl.append(xi == symbol(xs.str()));
            czRepl.append(xi == symbol(zs.str()));
            count++;
        }
        
        auto pls = get_pl_from(expr);
        int npls = pls.size()>0 ? ex_to<numeric>(pls[pls.size()-1].subs(lst{PL(wild())==wild()})).to_int() : -1;
        for(int i=0; i<npls+1; i++) {
            ostringstream pl;
            pl << "pl[" << i << "]";
            plRepl.append(PL(i) == symbol(pl.str()));
        }
        
        exset ftset;
        expr.find(FT(wild()), ftset);
        lst ftlst;
        for(auto it=ftset.begin(); it!=ftset.end(); it++) ftlst.append(*it);
        expr = expr.collect(ftlst);
        vector<pair<ex,ex>> ft_expr;
        ex ft;
        vector<ex> ftx;
        for(auto item : ftlst) {
            auto ftmp = item.subs(FT(wild())==wild());
            ft_expr.push_back(make_pair(ftmp, expr.coeff(item)));
            auto xys = get_xy_from(ftmp);
            if(xys.size()>=ftx.size()) {
                ft = ftmp;
                ftx = xys;
            }
        }
        
        bool need_contour_deformation = false;
        if(ft.has(x(wild()))) {
            auto tmp = ft.subs(nReplacements).expand();
            if(is_a<add>(tmp)) {
                need_contour_deformation = false;
                for(auto item : tmp) {
                    if(item.subs(x(wild())==1) < 0) {
                        need_contour_deformation = true;
                        break;
                    }
                }
                if(!need_contour_deformation) ft = 1; //note the difference with SDPrepare
            } else {
                if(tmp.subs(x(wild())==1) < 0) need_contour_deformation = true;
            }
        }
        
        ostringstream cppfn;
        cppfn << pid << "/" << rid << ".cpp";
        std::ofstream ofs;
        ofs.open(cppfn.str(), ios::out);
        if (!ofs) throw runtime_error("failed to open *.cpp file!");
        
/****************************************************************/
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
dCOMPLEX MatDetD(dCOMPLEX mat[], int n);
qCOMPLEX MatDetQ(qCOMPLEX mat[], int n);

#define expt(a,b) pow(a,b)
#define recip(a) pow(a,-1)

inline qREAL pow(qREAL x, qREAL y) { return powq(x, y); }
inline qREAL log(qREAL x) { return logq(x); }
inline qCOMPLEX pow(qCOMPLEX x, qREAL y) { return cpowq(x, y); }
inline qCOMPLEX log(qCOMPLEX x) { return clogq(x); }

#define Pi 3.1415926535897932384626433832795028841971693993751L
#define Euler 0.57721566490153286060651209008240243104215933593992L
#define iEpsilon complex<long double>(0,1.E-15L)

)EOF" << endl;
/****************************************************************/

        // long double
        auto cppL =  CppFormat(ofs, "L");
        ofs << "extern \"C\" " << endl;
        ofs << "int SD_D"<<rid<<"(const unsigned int xn, const qREAL xx[], const unsigned int yn, qREAL y[], const qREAL qpl[], const qREAL la[]) {" << endl;
        ofs << "dREAL x[xn];" << endl;
        ofs << "for(int i=0; i<xn; i++) x[i] = xx[i];" << endl;
        ofs << "dREAL pl["<<(npls<0 ? 1 : npls+1)<<"];" << endl;
        ofs << "for(int i=0; i<"<<(npls<0 ? 1 : npls+1)<<"; i++) pl[i] = qpl[i];" << endl;
        if(need_contour_deformation) {
            vector<symbol> ilas;
            for(int ii=0; ii<ftx.size(); ii++) {
                ostringstream ilaos;
                ilaos << "ila[" << ii << "]";
                ilas.push_back(symbol(ilaos.str()));
            }
            ofs << "dCOMPLEX z[xn];" << endl;
            ofs << "dCOMPLEX yy=0;" << endl;
            ofs << "dCOMPLEX ytmp, det;" << endl;
            ofs << "dCOMPLEX ila["<<ftx.size()<<"];" << endl;
            ofs << "dCOMPLEX mat["<<(ftx.size()*ftx.size())<<"];" << endl;
            ofs << "for(int i=0; i<"<<ftx.size()<<"; i++) ila[i] = complex<long double>(0, la[i]);" << endl;

            symbol s;
            for(auto &kv : ft_expr) {
                ofs << "for(int i=0; i<xn; i++) z[i] = x[i];" << endl;
                lst x0Repl;
                for(int ii=0; ii<ftx.size(); ii++) {
                    if(!kv.first.has(ftx[ii])) x0Repl.append(ftx[ii]==0);
                }
                
                for(int ii=0; ii<ftx.size(); ii++) {
                    auto xi = ftx[ii];
                    ofs << xi.subs(czRepl) << " = ";
                    auto zi = xi-xi*(1-xi)*ilas[ii]*(ft.subs(xi==s).diff(s).subs(s==xi).subs(x0Repl));
                    zi.subs(cxRepl).subs(plRepl).print(cppL);
                    ofs << ";" << endl;
                    for(int jj=0; jj<ftx.size(); jj++) {
                        ofs << "mat["<<(ii*ftx.size()+jj)<<"] = ";
                        zi.subs(ftx[jj]==s).diff(s).subs(s==ftx[jj]).subs(cxRepl).subs(plRepl).print(cppL);
                        ofs << ";" << endl;
                    }
                }
                ofs  << "det = MatDetD(mat, "<<ftx.size()<<");" << endl;
                ofs << "ytmp = ";
                kv.second.subs(czRepl).subs(plRepl).print(cppL);
                ofs << ";" << endl;
                ofs << "yy += det * ytmp;" << endl;
            }
        } else {
            auto tmp = expr.subs(FT(wild())==1).subs(cxRepl).subs(plRepl);
            ofs << "dCOMPLEX yy = ";
            tmp.print(cppL);
            ofs << ";" << endl;
        }
        ofs << "y[0] = yy.real();" << endl;
        ofs << "y[1] = yy.imag();" << endl;
        ofs << "return 0;" << endl;
        ofs << "}" << endl;
/****************************************************************/
ofs << R"EOF(
#undef Pi
#undef Euler
#undef iEpsilon
#define Pi 3.1415926535897932384626433832795028841971693993751Q
#define Euler 0.57721566490153286060651209008240243104215933593992Q
#define iEpsilon 1.E-30Qi
)EOF" << endl;
/****************************************************************/
        // Quadruple
        auto cppQ = CppFormat(ofs, "Q");
        ofs << "extern \"C\" " << endl;
        ofs << "int SD_Q"<<rid<<"(const unsigned int xn, const qREAL x[], const int unsigned yn, qREAL y[], const qREAL pl[], const qREAL la[]) {" << endl;
        if(need_contour_deformation) {
            vector<symbol> ilas;
            for(int ii=0; ii<ftx.size(); ii++) {
                ostringstream ilaos;
                ilaos << "ila[" << ii << "]";
                ilas.push_back(symbol(ilaos.str()));
            }
            ofs << "qCOMPLEX z[xn];" << endl;
            ofs << "qCOMPLEX yy=0;" << endl;
            ofs << "qCOMPLEX ytmp, det;" << endl;
            ofs << "qCOMPLEX ila["<<ftx.size()<<"];" << endl;
            ofs << "qCOMPLEX mat["<<(ftx.size()*ftx.size())<<"];" << endl;
            ofs << "for(int i=0; i<"<<ftx.size()<<"; i++) ila[i] = la[i] * 1.Qi;" << endl;

            symbol s;
            for(auto &kv : ft_expr) {
                ofs << "for(int i=0; i<xn; i++) z[i] = x[i];" << endl;
                lst x0Repl;
                for(int ii=0; ii<ftx.size(); ii++) {
                    if(!kv.first.has(ftx[ii])) x0Repl.append(ftx[ii]==0);
                }
                for(int ii=0; ii<ftx.size(); ii++) {
                    auto xi = ftx[ii];
                    ofs << xi.subs(czRepl) << " = ";
                    auto zi = xi-xi*(1-xi)*ilas[ii]*(ft.subs(xi==s).diff(s).subs(s==xi).subs(x0Repl));
                    zi.subs(cxRepl).subs(plRepl).print(cppQ);
                    ofs << ";" << endl;
                    for(int jj=0; jj<ftx.size(); jj++) {
                        ofs << "mat["<<(ii*ftx.size()+jj)<<"] = ";
                        zi.subs(ftx[jj]==s).diff(s).subs(s==ftx[jj]).subs(cxRepl).subs(plRepl).print(cppQ);
                        ofs << ";" << endl;
                    }
                }
                ofs  << "det = MatDetQ(mat, "<<ftx.size()<<");" << endl;
                ofs << "ytmp = ";
                kv.second.subs(czRepl).subs(plRepl).print(cppQ);
                ofs << ";" << endl;
                ofs << "yy += det * ytmp;" << endl;
            }
        } else {
            auto tmp = expr.subs(FT(wild())==1).subs(cxRepl).subs(plRepl);
            ofs << "qCOMPLEX yy = ";
            tmp.print(cppQ);
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
        remove(cppfn.str().c_str());
        return lst{ rid, xs.size(), kv.first, ft };
    }, "int", Verbose, false);
    
    ostringstream sofn, garfn, cmd;
    if(key != NULL) {
        sofn << key << ".so";
        garfn << key << ".pref.gar";
        lst gar_res;
        for(auto &item : res) gar_res.append(item);
        archive ar;
        ar.archive_ex(gar_res, "res");
        ar.archive_ex(19790923, "c");
        ofstream out(garfn.str());
        out << ar;
        out.close();
    } else {
        sofn << pid << ".so";
        for(auto &item : res) soIntegrands.push_back(ex_to<lst>(item));
    }
    
    CompileMatDet();
    cmd << "g++ -fPIC -shared -lquadmath " << CFLAGS << " -o " << sofn.str() << " " << pid << "/*.o";
    system(cmd.str().c_str());
    cmd.clear();
    cmd.str("");
    cmd << "rm -rf " << pid;
    system(cmd.str().c_str());
}

// need Parameter
void SD::ContourPrepares(const char * key) {
    if(IsZero) return;
    
    if(key != NULL) {
        ostringstream garfn;
        garfn << key << ".pref.gar";
        archive ar;
        ifstream in(garfn.str());
        in >> ar;
        in.close();
        auto c = ar.unarchive_ex(ParallelSymbols, "c");
        auto res = ar.unarchive_ex(ParallelSymbols, "res");
        if(c!=19790923) throw runtime_error("*.pref.gar error!");
        for(auto item : ex_to<lst>(res)) soIntegrands.push_back(ex_to<lst>(item));
    }
    
    lst fts;
    bool noFT = true;
    for(auto item : soIntegrands) {
        if(item.has(x(wild()))) {
            fts.append(item.op(3));
            noFT = false;
        }
    }
    fts.sort();
    fts.unique();
    
    lst plRepl;
    for(auto kv : Parameter) plRepl.append(PL(kv.first)==kv.second);
    fts=ex_to<lst>(subs(fts, plRepl));
    if(fts.has(PL(wild()))) {
        cout << "fts has PL: " << fts << " ; " << plRepl  << endl;
        assert(false);
    }
    
    if(noFT) return;
    
    if(Verbose > 0) cout << now() << " - ContourPrepares ..." << endl << flush;
    
    vector<ex> ftvec;
    for(auto item : fts) ftvec.push_back(item);
    
    auto ppid = getpid();
    ostringstream cmd;
    cmd << "mkdir -p " << ppid;
    system(cmd.str().c_str());
    
    vector<ex> res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, ftvec, [&](auto &ft, auto rid) {
        // return lst{ F-term, lst{lambda-i, lambda-max} }
        // with I*[lambda-i]*lambda, lambda < lambda-max
        auto xs = get_xy_from(ft);
        lst las, dfs;
        if(xs.size()<1) {
            las.append(numeric(100));
            return lst{ ft, las };
        }
        
        symbol s;
        double dfmax = -1;
        for(int i=0; i<xs.size(); i++) {
            auto df = ft.subs(xs[i]==s).diff(s).subs(s==xs[i]);
            dfs.append(df);
            auto max_df = -FindMinimum(-fabs(df));
            las.append(max_df);
            if(dfmax < max_df) dfmax = max_df;
        }
        
        ex zs[xs.size()];
        symbol ila("ila");
        for(int i=0; i<xs.size(); i++) {
            //TODO: add other schema
            auto tmp = 1/dfmax;
            auto dla = ex_to<numeric>(las.op(i)).to_double();
            if(dla > 1E-2 * dfmax) tmp = 1/dla;
            
            las.let_op(i) = tmp;
            zs[i] = xs[i] - xs[i]*(1-xs[i])*dfs.op(i)*ila*tmp;
        }
        
        ostringstream cppfn, sofn;
        auto pid = getpid();
        cppfn << ppid << "/" << pid << "-la.cpp";
        sofn << ppid << "/" << pid << "-la.so";
        std::ofstream ofs;
        ofs.open(cppfn.str(), ios::out);
        if (!ofs) throw runtime_error("failed to open *.cpp file!");
        
        lst cxRepl, czRepl;
        int count = 0;
        for (auto xi : xs) {
            ostringstream xs, zs;
            xs << "x[" << count << "]";
            cxRepl.append(xi == symbol(xs.str()));
            zs << "z[" << count << "]";
            czRepl.append(xi == symbol(zs.str()));
            count++;
        }

/****************************************************************/
ofs << R"EOF(
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <iostream>
using namespace std;

#define Pi 3.1415926535897932384626433832795028841971693993751L
#define Euler 0.57721566490153286060651209008240243104215933593992L

typedef long double dREAL;
typedef complex<long double> dCOMPLEX;

#define expt(a,b) pow(a,b)
#define recip(a) pow(a,-1)
)EOF" << endl;
/****************************************************************/
        auto cppL = CppFormat(ofs, "L");
        ofs << "extern \"C\" " << endl;
        ofs << "dREAL minLaFunc(int xn, dREAL* x) {" << endl;
        ofs << "dCOMPLEX ila = complex<long double>(0, x[xn-1]);" <<endl;
        ofs << "dCOMPLEX z[xn];" << endl;
        for(int i=0; i<xs.size(); i++) {
            ofs << "z["<<i<<"] = ";
            zs[i].subs(cxRepl).print(cppL);
            ofs << ";" << endl;
        }
        ofs << "dCOMPLEX zf = " << endl;
        ft.subs(czRepl).print(cppL);
        ofs << ";" << endl;
        ofs << "return -zf.imag();" << endl; // find max image part, check with 0
        ofs << "}" << endl;
        
        cmd.clear();
        cmd.str("");
        cmd << "g++ -fPIC -shared " << CFLAGS << " -o " << sofn.str() << " " << cppfn.str();
        system(cmd.str().c_str());
        
        void* module = nullptr;
        module = dlopen(sofn.str().c_str(), RTLD_NOW);
        if (module == nullptr) throw std::runtime_error("could not open compiled module!");
        auto fp = (MinimizeBase::FunctionType)dlsym(module, "minLaFunc");
        assert(fp!=NULL);
        
        dREAL laBegin = 0, laEnd = 100, min;
        int nvars = xs.size()+1;
        dREAL UB[nvars];
        for(int i=0; i<nvars; i++) UB[i] = 1;
        
        min = laEnd;
        while(true) {
            UB[nvars-1] = min;
            dREAL res = Minimizer->FindMinimum(nvars, fp, UB);
            
            if(res < -5E-4) laEnd = min;
            else laBegin = min;
            
            if(laEnd - laBegin < 1E-3 * laEnd) break;
            min = (laBegin + laEnd) / 2.0;
        }
        min = laBegin;
        
        if(use_dlclose) dlclose(module);
        cmd.clear();
        cmd.str("");
        cmd << "rm " << cppfn.str() << " " << sofn.str();
        system(cmd.str().c_str());
                
        las.append(numeric((double)min));
        return lst{ ft, las };
    
    }, "la", Verbose);
    
    ostringstream garfn;
    if(key != NULL) {
        garfn << key << ".la.gar";
        lst gar_res;
        for(auto &item : res) gar_res.append(item);
        archive ar;
        ar.archive_ex(gar_res, "res");
        ar.archive_ex(19790923, "c");
        ofstream out(garfn.str());
        out << ar;
        out.close();
    } else {
        for(auto &item : res) LambdaFT[item.op(0)] = item.op(1);
    }
    
    cmd.clear();
    cmd.str("");
    cmd << "rm -rf " << ppid;
    system(cmd.str().c_str());
}

// need Parameter
void SD::Integrates(const char * key) {
    if(IsZero) return;
    
    if(Verbose > 0) cout << now() << " - Integrates ..." << endl << flush;
    
    auto pid = getpid();
    ostringstream sofn, cmd;
    if(key == NULL) {
        sofn << pid << ".so";
    } else {
        sofn << key << ".so";
        ostringstream garfn;
        garfn << key << ".pref.gar";
        archive ar;
        ifstream in(garfn.str());
        in >> ar;
        in.close();
        auto c = ar.unarchive_ex(ParallelSymbols, "c");
        auto res = ar.unarchive_ex(ParallelSymbols, "res");
        if(c!=19790923) throw runtime_error("*.pref.gar error!");
        for(auto item : ex_to<lst>(res)) soIntegrands.push_back(ex_to<lst>(item));
        
        garfn.clear();
        garfn.str("");
        garfn << key << ".la.gar";
        if(file_exists(garfn.str().c_str())) {
            archive la_ar;
            ifstream la_in(garfn.str());
            la_in >> la_ar;
            la_in.close();
            auto la_c = la_ar.unarchive_ex(ParallelSymbols, "c");
            auto la_res = la_ar.unarchive_ex(ParallelSymbols, "res");
            if(la_c!=19790923) throw runtime_error("*.pref.gar error!");
            for(auto item : ex_to<lst>(la_res)) {
                LambdaFT[item.op(0)] = item.op(1);
            }
        }
    }
    
    void* module = nullptr;
    module = dlopen(sofn.str().c_str(), RTLD_NOW);
    if (module == nullptr) throw std::runtime_error("could not open compiled module!");
    if(key == NULL) remove(sofn.str().c_str());
    
    int nprar = 0;
    lst plRepl;
    for(auto kv : Parameter) {
        plRepl.append(PL(kv.first)==kv.second);
        if(kv.first>nprar) nprar = kv.first;
    }
    
    int total = soIntegrands.size(), current = 0;
    ResultError = 0;
    for(auto &item : soIntegrands) {
        if(Verbose > 1) {
            cout << "\r                                      \r" << flush;
            cout << "\r  \\--Evaluating [" <<(++current)<<"/"<<total<< "] ... " << flush;
        }
        
        unsigned int xsize = ex_to<numeric>(item.op(1)).to_int();
        if(xsize<1) {
            Digits = 35;
            ResultError +=  VE(item.op(0).subs(plRepl).evalf(),0) * item.op(2).subs(plRepl);
            continue;
        }
        
        int rid = ex_to<numeric>(item.op(0)).to_int();
        auto co = item.op(2).subs(plRepl).expand().collect(lst{eps, ep}, true);
        if(co.is_zero()) continue;
        assert(!co.has(PL(wild())));
        qREAL cmax = -1;
        for(int si=co.ldegree(eps); si<=co.degree(eps); si++) {
            auto tmp = co.coeff(eps, si);
            assert(!tmp.has(eps));
            for(int i=tmp.ldegree(ep); i<=tmp.degree(ep); i++) {
                auto ccRes = tmp.coeff(ep, i).evalf().expand();
                lst css;
                css.append(ccRes);
                if(is_a<add>(ccRes)) {
                    for(auto item : ccRes) css.append(item);
                }
                
                for(int ci=0; ci<css.nops(); ci++) {
                    auto nt = css.op(ci).subs(nReplacements).evalf();
                    nt = abs(nt).evalf(); // no PL here, then nReplacements
                    if(!is_a<numeric>(nt)) {
                        cout << "nt: " << nt << endl;
                        assert(false);
                    }
                    assert(!nt.has(ep));
                    qREAL qnt = CppFormat::ex2q(nt);
                    if(qnt > cmax) cmax = qnt;
                }
            }
        }
        if(cmax<=0) {
            cout << "cmax<=0 with co = " << co <<endl;
            assert(false);
        }
        
        ostringstream fname;
        fname << "SD_D" << rid;
        auto fp = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
        assert(fp!=NULL);
        fname.clear();
        fname.str("");
        fname << "SD_Q" << rid;
        auto fpQ = (IntegratorBase::SD_Type)dlsym(module, fname.str().c_str());
        assert(fpQ!=NULL);
        
        auto las = LambdaFT[item.op(3).subs(plRepl)];
        qREAL lambda[las.nops()];
        qREAL paras[nprar+1];
        for(auto kv : Parameter) paras[kv.first] = CppFormat::ex2q(kv.second);
        
        
        if(is_a<lst>(las)) {
            qREAL lamax = CppFormat::ex2q(las.op(las.nops()-1));
            if(lamax > LambdaMax) lamax = LambdaMax;
            
            Integrator->RunMAX = 1;
            Integrator->RunPTS = TryPTS;
            Integrator->EpsAbs = 0;
            Integrator->EpsRel = 0;

            int smin = -1;
            for(int stry=0; stry<5; stry++) {
                smin = -1;
                ex emin = 0;
                for(int s=0; s<LambdaSplit; s++) {
                    for(int i=0; i<las.nops()-1; i++) {
                        lambda[i] = (s+1) * lamax *  CppFormat::ex2q(las.op(i)) / LambdaSplit;
                    }
                    auto res = Integrator->Integrate(xsize, fp, fpQ, paras, lambda);
                    if(res.has(NaN)) continue;
                    auto res_tmp = res.subs(VE(wild(1), wild(2))==wild(2));
                    auto err = real_part(res_tmp);
                    if(err < imag_part(res_tmp)) err = imag_part(res_tmp);
                    if(smin<0 || err < emin) {
                        smin = s;
                        emin = err;
                        if(emin < CppFormat::q2ex(EpsAbs)) {
                            smin = -2;
                            ResultError += co * res;
                            break;
                        }
                    }
                } 
                if(smin == -2) break; 
                
                if(smin == -1) {
                    cout << "smin = -1, optimized lambda NOT found!" << endl;
                    assert(false);
                }
                
                if(smin < 2) {
                    lamax = (smin+1) * lamax / LambdaSplit;
                    continue;
                } else break;
            }
            
            if(smin == -2) continue;
            
            for(int i=0; i<las.nops()-1; i++) {
                lambda[i] = (smin+1) * lamax * CppFormat::ex2q(las.op(i)) / LambdaSplit;
            }
        }
        
        Integrator->RunMAX = RunMAX;
        Integrator->RunPTS = RunPTS;
        Integrator->EpsAbs = EpsAbs/cmax;
        Integrator->EpsRel = 0;
        auto res = Integrator->Integrate(xsize, fp, fpQ, paras, lambda);
        if(res.has(NaN)) ResultError += NaN;
        else ResultError += co * res;
    }
    
    if(use_dlclose) dlclose(module);
    if(total>0 && Verbose > 1) cout << "@" << now(false) << endl;
    
    ResultError = VESimplify(ResultError, epN);
}

void SD::Initialize(XIntegrand xint) {
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

void SD::Evaluate(FeynmanParameter fp) {    
    cout << endl << "Starting @ " << now() << endl;
    if(SecDec==NULL) SecDec = new SecDecG();
    if(Integrator==NULL) Integrator = new HCubature();
    
    //if(Minimizer==NULL) Minimizer = new HookeJeeves();
    if(Minimizer==NULL) Minimizer = new MinUit();
    
    Initialize(fp);
    if(FunExp.size()<1) return;
    Scalelesses();
        
    RemoveDeltas();
    SDPrepares();
    EpsEpExpands();
    IntPrepares();
    ContourPrepares();
    Integrates();
    delete SecDec;
    delete Integrator;
    delete Minimizer;
    cout << "Finished @ " << now() << endl << endl;
}

void SD::Evaluate(XIntegrand xint) {
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
    IntPrepares();
    ContourPrepares();
    Integrates();
    delete SecDec;
    delete Integrator;
    delete Minimizer;
    cout << "Finished @ " << now() << endl << endl;
}

double SD::FindMinimum(ex expr) {
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
    
/****************************************************************/
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

#define expt(a,b) pow(a,b)
#define recip(a) pow(a,-1)
)EOF" << endl;
/****************************************************************/

    auto cppL = CppFormat(ofs, "L");
    ofs << "extern \"C\" " << endl;
    ofs << "dREAL minFunc(int xn, dREAL* x) {" << endl;
    auto tmp = expr.subs(cxRepl);
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
    
    double min = Minimizer->FindMinimum(count, fp, UB, LB);
    
    if(use_dlclose) dlclose(module);
    cmd.clear();
    cmd.str("");
    cmd << "rm " << cppfn.str() << " " << sofn.str();
    system(cmd.str().c_str());
    
    return min;
}

int SD::epRank(ex expr) {
    if(!expr.has(ep)) return 0;
    int p = -5;
    while(true) {
        auto tmp = normal(series_to_poly(expr.series(ep, p)));
        if(!tmp.is_zero()) {
            return tmp.ldegree(ep);
        } else p++;
    }
    throw runtime_error("epRank error!");
}

int SD::epsRank(ex expr) {
    if(!expr.has(eps)) return 0;
    int p = -5;
    while(true) {
        auto tmp = normal(series_to_poly(expr.series(eps, p)));
        if(!tmp.is_zero()) {
            return tmp.ldegree(eps);
        } else p++;
    }
    throw runtime_error("epsRank error!");
}

ex SD::VEResult() {
    return HepLib::VEResult(ResultError);
}

ex SD::PrefactorFIESTA(int nLoop) {
    return  pow(I*pow(Pi,2-ep)*exp(-ep*Euler), -ex(nLoop));
}

}
