#include "MB.h"
#include <math.h>
#include <cmath>

namespace HepLib {

symbol const MB::ep("ep");

void Replacement2(exmap &repl) {
    auto tmp = repl;
    for(auto &kv : repl) {
        kv.second = kv.second.subs(tmp);
    }
}

void MB::Initialize(FeynmanParameter fp) {

    if(fp.Propagator.nops() != fp.Exponent.nops()) {
        cerr << ErrColor << "the length of Propagator and Exponent are NOT equal." << RESET << endl;
        exit(1);
    }
            
    auto ps = fp.Propagator;
    auto ns = fp.Exponent;
    auto ls = fp.LoopMomenta;
    
    Replacement2(fp.lReplacement);
    Replacement2(fp.nReplacement);
    
    auto lsubs = fp.lReplacement;
    auto nsubs = fp.nReplacement;
    auto nReplacement = fp.nReplacement;
    
    if(Verbose > 1) cout << "  Initialize @ " << now() << endl;
    
    ex asgn = 1;
    ex a = 0;
    ex ad = (4-2*ep)*ls.nops();
    int xn = ps.nops();
    ex rem = 0;
    exmap xtNeg;
    
    ex pre = fp.Prefactor; // come from below
    for(int i=0; i<ps.nops(); i++) {
        bool ltQ = false; {
            auto tps = ps.op(i).expand().subs(lsubs);
            for(auto lsi : ls) {
                if(tps.has(lsi)) {
                    ltQ = true;
                    break;
                }
            }
        }
        
        if(ltQ) a += ns.op(i);
        ex sgn = 0;
        
        if(!ltQ) {
            pre = pre * pow(ps.op(i).expand().subs(lsubs), ex(0)-ns.op(i));
            ns[i] = 0;
            ps[i] = 1;
            continue;
        } else if(ns.op(i).info(info_flags::negint) && is_zero(ns.op(i))) {
            xtNeg[x(i)]=0;
            if(!is_zero(ns.op(i))) {
                asgn *= pow(-1, ns.op(i));
                rem += x(i) * ps.op(i);
            }
            continue;
        }

        auto p = ps.op(i).expand().subs(lsubs).subs(nsubs);
        p = p.subs(lsubs).subs(nsubs);
        // check loop^2
        for(auto m : ls) {
            if(!is_a<numeric>(p.coeff(m,2))) {
                cout << ErrColor << "not numeric: " << p.coeff(m,2) << endl;
                cout << "nsubs = " << nsubs << RESET << endl;
                exit(1);
            }
            numeric nm = ex_to<numeric>(p.coeff(m,2));
            if(nm.is_zero()) continue;
            sgn = nm>0 ? -1 : 1;
            break;
        }
        
        // others
        if(sgn.is_zero()) {
            sgn = 1;
            if(is_a<numeric>(p) && ex_to<numeric>(p)>0) sgn = -1;
            cout << " - Warning: Can NOT determine the iEpsilon sign." << endl;
            cout << " - " << p << " from " << ps.op(i) << endl;
        }
        
        p = (ps.op(i)*sgn);
        if(sgn==-1) asgn *= exp(I * Pi * ns.op(i)); // sgn<0: +iep, sgn>0: -iep, so I*Pi here
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
        rem = normal(rem.subs(lsubs).subs(lsubs));
        u = normal(u.subs(lsubs));
        
        cu *= u;
        auto u_nd = numer_denom(u);
        ex usgn = u_nd.op(1).subs(xtNeg).subs(x(w)==ex(1)/2).subs(nsubs);
        if(usgn.is_zero()) usgn = u_nd.op(1).subs(xtNeg).subs(x(w)==ex(1)/3).subs(nsubs);
        if(usgn.is_zero()) {
            cerr << ErrColor << "usgn is zero!" << RESET << endl;
            exit(1);
        }
        usgn = normal(usgn)>0 ? 1 : -1;
        
        if(!xPositive(normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacement).subs(lst{
            CV(w1,w2)==w2, ep==ex(1)/111
        }))) {
            cerr <<ErrColor << "NOT positive - un: " << normal(usgn*u_nd.op(0)).subs(xtNeg).subs(nReplacement).subs(lst{
                CV(w1,w2)==w2, ep==ex(1)/111
            }) << RESET << endl;
            exit(1);
        }
        if(!xPositive(normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacement).subs(lst{
            CV(w1,w2)==w2, ep==ex(1)/111
        }))) {
            cerr << ErrColor << "NOT positive - ud: " << normal(usgn*u_nd.op(1)).subs(xtNeg).subs(nReplacement).subs(lst{
                CV(w1,w2)==w2, ep==ex(1)/111
            }) << RESET << endl;
            exit(1);
        }
        
        uList1.append(usgn*u_nd.op(0));
        uList2.append(-(4-2*ep)/2);
        if((usgn*u_nd.op(0)) != 1) {
            uList1.append(usgn*u_nd.op(1));
            uList2.append((4-2*ep)/2);
        }
    } else {
        rem = normal(rem.subs(lsubs).subs(lsubs));
    }

    u = normal(cu);
    auto u_nd = numer_denom(u);
    rem = normal(rem * u);
    auto rem_nd = numer_denom(rem);
    
    ex usgn = u_nd.op(1).subs(xtNeg).subs(x(w)==ex(1)/2).subs(nsubs);
    usgn = normal(usgn)>0 ? 1 : -1;
    ex fsgn = rem_nd.op(1).subs(xtNeg).subs(x(w)==ex(1)/2).subs(nsubs);
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
    if(ns.op(i).info(info_flags::negint)) {
        for(int j=0; j<-ns.op(i); j++) {
            vector<pair<lst, lst>> nret;
            for(auto kv : ret) {
                auto plst = kv.first;
                auto nlst = kv.second;
                for(int ij=0; ij<nlst.nops(); ij++) {
                    auto dtmp = nlst.op(ij) * diff_ex(plst.op(ij),x(i));
                    if(dtmp.is_zero()) continue;
                    auto plst2 = plst;
                    auto nlst2 = nlst;
                    if((nlst.op(ij)-1).is_zero()) {
                        plst2[ij] = dtmp;
                    } else {
                        nlst2[ij] = nlst.op(ij)-1;
                        int nn = plst.nops();
                        if(!(nlst.op(nn-1)-1).is_zero()) {
                            plst2.append(dtmp);
                            nlst2.append(1);
                        } else plst2[nn-1] = plst.op(nn-1) * dtmp;
                    }
                    nret.push_back(make_pair(plst2, nlst2));
                }
            }
            ret = nret;
        }
        
        for(auto &kv : ret) {
            lstHelper::map_inplace(kv.first, [](auto &&e) { return e.subs(x(i)==0); });
            lstHelper::map_inplace(kv.second, [](auto &&e) { return e.subs(x(i)==0); });
        }
    }}

    // simplification
    // ex pre = fp.Prefactor; // moved to above
    pre *= asgn * pow(I,ls.nops()) * pow(Pi, ad/2) * tgamma(a-ad/2);
    for(int i=0; i<ns.nops(); i++) {
        if(is_a<numeric>(ns.op(i)) && ns.op(i)<=0) continue;
        pre /= tgamma(ns.op(i));
    }
    
    ex xpre = 1;
    for(int i=0; i<ns.nops(); i++) {
        if(is_a<numeric>(ns.op(i)) && ns.op(i)<=1) continue;
        else {
            for(auto &kv : ret) {
                kv.first.append(x(i));
                kv.second.append(ns.op(i)-1);
             }
        }
    }

    for(auto &kv : ret) {
        kv.first.append(pre);
        kv.second.append(1);
        if(xpre != 1) {
            kv.first.append(xpre);
            kv.second.append(1);
        }
        lstHelper::map_inplace(kv.first, [](auto &&e) { return collect_common_factors(e); });
        lstHelper::map_inplace(kv.second, [](auto &&e) { return collect_common_factors(e); });
    }

    lst delta;
    for(int i=0; i<ns.nops(); i++) {
        if(is_a<numeric>(ns.op(i)) && ns.op(i)<=0) continue;
        delta.append(x(i));
    }
    Deltas.push_back(delta);
    FunExp = ret;
    
}



}
