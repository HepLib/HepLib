#include "SD.h"
#include <math.h>
#include <cmath>

namespace HepLib {

symbol const SD::ep("ep");
symbol const SD::eps("eps");
symbol const SD::vs("s");
symbol const SD::vz("z");
symbol const SD::iEpsilon("iEpsilon");
realsymbol const SD::NaN("NaN");
bool SD::use_dlclose = true;
bool SD::debug = false;

vector<exmap> SecDecBase::x2y(const lst &xpols, bool all_in_one) {
    if(all_in_one) {
        ex xpol = 1;
        for(auto item : xpols) xpol *= item;
        return x2y(xpol);
    }
    
    vector<exmap> vec_map;
    exmap map0;
    map0[x(-1)] = 1;
    auto xs = get_x_from(xpols);
    for(auto xi : xs) map0[xi] = xi.subs(x(wild())==y(wild()));
    vec_map.push_back(map0);
    
    for(int i=0; i<xpols.nops(); i++) {
        auto xpol = xpols.op(i);
        vector<exmap> vec_map2;
        for(auto map : vec_map) {
            auto cxpol = xpol.subs(map).subs(y(wild())==x(wild()));
            auto vec_map3 = x2y(cxpol);
            for(auto map3 : vec_map3) {
                exmap new_map;
                lst xy_lst;
                for(auto kv : map) {
                    auto tmp = kv.second.subs(y(wild())==x(wild())).subs(map3);
                    if(kv.first==x(-1)) {
                        tmp *= map3[x(-1)];
                    }
                    new_map[kv.first] = tmp;
                    xy_lst.append(tmp);
                }
                
                // handle remaining x's
                if(xy_lst.has(x(wild()))) {
                    vector<int> y_free_index;
                    for(int j=0; j<xs.size()+10; j++) {
                        if(!xy_lst.has(y(j))) y_free_index.push_back(j);
                    }
                    auto xs = get_x_from(xy_lst);
                    lst xRepl;
                    assert(xs.size()<=y_free_index.size());
                    for(int j=0; j<xs.size(); j++) {
                        xRepl.append(xs[j]==y(y_free_index[j]));
                    }
                    for(auto kv : new_map) new_map[kv.first] = kv.second.subs(xRepl);
                }
                
                vec_map2.push_back(new_map);
            }
        }
        vec_map = vec_map2;
    }
    return vec_map;
}

/*-----------------------------------------------------*/
/*                   CheckAtEnd                        */
/*-----------------------------------------------------*/
bool SD::IsBad(ex f, vector<exmap> vmap) {
    for(auto &vi : vmap) {
        auto ft = f.subs(vi);
        auto xs_tmp = get_x_from(ft);
        auto ys_tmp = get_y_from(ft);
        int ysn = ys_tmp.size();
        vector<int> y_free_index;
        for(int j=0; j<ys_tmp.size()+10; j++) {
            if(!ft.has(y(j))) y_free_index.push_back(j);
        }
        assert(xs_tmp.size()<=y_free_index.size());
        for(int i=0; i<xs_tmp.size(); i++) {
            vi[xs_tmp[i]] = y(y_free_index[i]);
        }
        if(xs_tmp.size()>0) ft = ft.subs(vi);

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

vector<lst> SD::AutoEnd(lst po_ex) {
    assert(po_ex.nops()<3);
    lst const exlist = ex_to<lst>(po_ex.op(1));
    assert((exlist.op(0)-1).is_zero());
    auto xs = get_x_from(po_ex.op(0));
    if(xs.size()<1) xs = get_y_from(po_ex.op(0));
    int nx = xs.size();
    for(int nn=0; nn<=nx; nn++) {
    for(int pi=0; pi<std::pow(2, nx); pi++) {
        int cpi = pi, cn1 = 0;
        for(int i=0; i<nx; i++) {
            if((cpi % 2) == 1) cn1++;
            cpi /= 2;
        }
        if(cn1 != nn) continue;
        
        lst polists = lst{ po_ex.op(0) };
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
            
            for(int ni=0; ni<polist.nops(); ni++) {
                auto po = polist.op(ni);
                auto expo = exlist.op(ni); // note that expo may decrease when using diff
                if(is_a<numeric>(expo) && ex_to<numeric>(expo).is_pos_integer()) continue;
                if(IsBad(po, vmap)) {
                    OK = false;
                    break;
                }
            }
            if(!OK) break;
        }

        if(OK) {
            vector<lst> res;
            for(auto item : polists) res.push_back(lst{ex_to<lst>(item), exlist});
            return res;
        }
    }}
    
    cerr << RED << "polynormial list: " << po_ex.op(0) << RESET << endl;
    cerr << RED << "AutoEnd Failed @ ALL possible bisections!" << RESET << endl;
    assert(false);
    return vector<lst>();
}

/*-----------------------------------------------------*/
/*                       SD                            */
/*-----------------------------------------------------*/
//return a lst, element pattern: { {{x1,n1}, {x2,n2}, ...}, {{e1, n1},{e2,n2}, ...} }
//e1 is a const term, e2 still the F-term
vector<lst> SD::DS(lst po_ex) {
    // 1st element in input polist is the constant term, guess NOT necessary
    // 2nd element in input polist is the F-term, required!
    lst const polist = ex_to<lst>(po_ex.op(0));
    lst const exlist = ex_to<lst>(po_ex.op(1));
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
        vector<int> y_free_index;
        for(int j=0; j<xs_tmp.size()+ys_tmp.size()+10; j++) {
            if(!ypolist.has(y(j))) y_free_index.push_back(j);
        }
        assert(xs_tmp.size()<=y_free_index.size());
        for(int i=0; i<xs_tmp.size(); i++) {
            vi[xs_tmp[i]] = y(y_free_index[i]);
        }
        if(xs_tmp.size()>0) ypolist = ypolist.subs(vi);
        assert(!ypolist.has(x(wild())));

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
            auto tmp = ft.subs(nReplacements).subs(lst{
                CV(wild(1),wild(2))==wild(2), ep==ex(1)/111, eps==ex(1)/1111
            }).expand();
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
            auto tmp = ft.subs(y(wild())==1).subs(nReplacements).subs(lst{
                CV(wild(1),wild(2))==wild(2), ep==ex(1)/111, eps==ex(1)/1111
            });
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
                        auto tr = item.subs(nReplacements).subs(lst{
                            CV(wild(1),wild(2))==wild(2), ep==ex(1)/111, eps==ex(1)/1111
                        });
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
                    cerr << endl << RED << k << "^(" << nv << ") Found, Other regularization needed!" << RESET << endl;
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

// 1st element in [output 1st] is the constant term
// 2nd element in both [input 1st] and [output 1st] is the F-term
lst SD::Normalize(const lst &input) {
    ex const_term = 1;
    lst plst, nlst;
    lst in_plst = ex_to<lst>(input.op(0));
    lst in_nlst = ex_to<lst>(input.op(1));
    for(int i=0; i<in_plst.nops(); i++) {
        if(i!=1 && (in_nlst.op(i).is_zero() || in_plst.op(i)==ex(1))) continue;
        if(i!=1 && !in_plst.op(i).has(x(wild())) && !in_plst.op(i).has(y(wild()))) {
            const_term *= pow(in_plst.op(i), in_nlst.op(i));
        } else {
            auto ptmp = in_plst.op(i);
            auto ntmp = in_nlst.op(i);
            if(is_exactly_a<mul>(ptmp)) {
                ex tmul = 1;
                for(int j=0; j<ptmp.nops(); j++) {
                    auto tmp = ptmp.op(j);
                    if(!tmp.has(x(wild())) && !tmp.has(y(wild()))) {
                        if(is_a<numeric>(ntmp) && ex_to<numeric>(ntmp).is_integer()) {
                            const_term *=  pow(tmp,ntmp);
                        } else if((tmp-vs).is_zero() || tmp.match(pow(vs,wild()))) {
                            const_term *=  pow(tmp,ntmp);
                        } else if(!tmp.has(PL(wild())) && !tmp.has(vs)) {
                            auto tr = tmp.subs(nReplacements).subs(lst{
                                CV(wild(1),wild(2))==wild(2), ep==ex(1)/111, eps==ex(1)/1111
                            });
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
                
                if(i==1) {
                    plst.prepend(tmul);
                    nlst.prepend(ntmp);
                } else if(tmul != 1) {
                    plst.append(tmul);
                    nlst.append(ntmp);
                }
            } else {
                if(i==1) {
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
    for(int i=0; i<nlst.nops(); i++) {
        if(i!=1) np[plst[i]] += nlst[i];
    }
    map<ex,int,ex_is_less> inp;
    for(int i=0; i<nlst.nops(); i++) {
        if(i==1) {
            plst_comb.append(plst[1]);
            nlst_comb.append(nlst[1]);
        } else {
            if(inp[plst[i]]!=0) continue;
            plst_comb.append(plst[i]);
            nlst_comb.append(np[plst[i]]);
            inp[plst[i]] = 1;
        }
    }
    
    return (input.nops()>2) ? lst{plst_comb, nlst_comb, input.op(2)} : lst{plst_comb, nlst_comb};
}

/*-----------------------------------------------------*/
/*               's Funtions in SD                     */
/*-----------------------------------------------------*/
// working with or without Deltas
void SD::XReOrders() {
    if(IsZero) return;
    if(Integrands.size()<1) {
        for(auto &fe : FunExp) {
            exset xset;
            find(fe, x(wild()), xset);
            
            vector<ex> xs;
            for(auto ii : xset) xs.push_back(ii);
            sort(xs.begin(), xs.end(), [&](const auto &a, const auto &b){
                return ex_to<numeric>(normal((b-a)).subs(lst{x(wild())==wild()})).is_positive();
            });
            
            lst x2y;
            for(int i=0; i<xs.size(); i++) {
                x2y.append(xs[i]==y(i));
            }
            
            fe = ex_to<lst>(subs(fe, x2y).subs(y(wild())==x(wild())));
        }
    } else {
        for(auto &vint : Integrands) {
            exset xset;
            find(vint, x(wild()), xset);
            
            vector<ex> xs;
            for(auto ii : xset) xs.push_back(ii);
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

// working with or without Deltas
void SD::Normalizes() {
    if(IsZero) return;

    vector<lst> funexp;
    for(auto fe : FunExp) {
        funexp.push_back(Normalize(fe));
    }
    FunExp.clear();
    FunExp.shrink_to_fit();
    
    exmap fn;
    for(auto fe : funexp) {
        ex key = 1;
        if(fe.nops()>2) key = WF(fe.op(2));
        for(int i=1; i<fe.op(0).nops(); i++) key *= pow(fe.op(0).op(i), fe.op(1).op(i));
        fn[key] += fe.op(0).op(0);
    }
    
    exmap ifn;
    for(auto fe : funexp) {
        ex key = 1;
        if(fe.nops()>2) key = WF(fe.op(2));
        for(int i=1; i<fe.op(0).nops(); i++) key *= pow(fe.op(0).op(i), fe.op(1).op(i));
        if(ifn[key]>0) continue;
        lst fun, exp;
        fun.append(fn[key]);
        exp.append(1);
        for(int i=1; i<fe.op(0).nops(); i++) {
            fun.append(fe.op(0).op(i));
            exp.append(fe.op(1).op(i));
        }
        if(fe.nops()>2) FunExp.push_back(lst{fun, exp, fe.op(2)});
        else FunExp.push_back(lst{fun, exp});
        ifn[key] = 1;
    }
}

// working with or without Deltas
void SD::XTogethers() {
    vector<lst> funexp;
    for(auto fe : FunExp) {
        funexp.push_back(fe);
    }
    FunExp.clear();
    FunExp.shrink_to_fit();
    
    map<ex, ex, ex_is_less> fe_cc;
    for(auto fe : funexp) {
        lst fun, exp;
        fun.append(1);
        exp.append(1);
        fun.append(fe.op(0).op(1));
        exp.append(fe.op(1).op(1));
        ex rem = 1;
        for(int i=0; i<fe.op(1).nops(); i++) {
            if(i==1) continue;
            auto expi = fe.op(1).op(i);
            if(is_a<numeric>(expi) && ex_to<numeric>(expi).is_nonneg_integer()) {
                rem *= pow(fe.op(0).op(i), fe.op(1).op(i));
            } else {
                fun.append(fe.op(0).op(i));
                exp.append(fe.op(1).op(i));
            }
        }
        auto key = fe;
        key.let_op(0) = fun;
        key.let_op(1) = exp;
        fe_cc[key] += rem;
    }
    
    for(auto kv : fe_cc) {
        lst fe = ex_to<lst>(kv.first);
        let_op_append(fe, 0, kv.second);
        let_op_append(fe, 1, 1);
        FunExp.push_back(fe);
    }
}

// working with or without Deltas
void SD::XExpands() {
    vector<lst> funexp;
    for(auto fe : FunExp) {
        funexp.push_back(fe);
    }
    FunExp.clear();
    FunExp.shrink_to_fit();
    
    for(auto fe : funexp) {
        lst fun = lst{fe.op(0).op(0), fe.op(0).op(1)};
        lst exp = lst{fe.op(1).op(0), fe.op(1).op(1)};
        ex rem = 1;
        for(int i=2; i<fe.op(1).nops(); i++) {
            auto expi = fe.op(1).op(i);
            if(is_a<numeric>(expi) && ex_to<numeric>(expi).is_nonneg_integer()) {
                rem *= pow(fe.op(0).op(i), fe.op(1).op(i));
            } else {
                fun.append(fe.op(0).op(i));
                exp.append(fe.op(1).op(i));
            }
        }
        
        rem = mma_collect(rem, x(wild()));
        
        lst rem_lst;
        if(is_a<add>(rem)) {
            for(auto item : rem) rem_lst.append(item);
        } else {
            rem_lst.append(rem);
        }
        
        for(auto item : rem_lst) {
            auto fe2 = fe;
            fe2.let_op(0) = fun;
            fe2.let_op(1) = exp;
            let_op_append(fe2, 0, item);
            let_op_append(fe2, 1, 1);
            FunExp.push_back(fe2);
        }
    }
    
}

// Section 2.1 @ https://arxiv.org/pdf/1712.04441.pdf
// also refers to Feng/Thinking.pdf
void SD::Scalelesses(bool verb) {
    if(IsZero) return;
    if(verb) cout << now() << " - Scaleless: " << FunExp.size() << " :> " << flush;

    vector<ex> sl_res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, FunExp, [&](auto &funexp, auto rid) {
        if(funexp.nops()<3) return funexp;
        symbol s;
        auto fun = funexp.op(0);
        auto exp = funexp.op(1);
        auto deltas = funexp.op(2);
        bool is0;
        for(int di=0; di<deltas.nops(); di++) {
            auto delta = ex_to<lst>(deltas.op(di));
            is0 = false;
            if(delta.nops()<2) continue;
            // to make sure the integrand is projective
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
                ex n_s = 0;
                for(int j=0; (j<delta.nops() && ci>0); j++) {
                    if((ci % 2)==1) {
                        sRepl.append(delta[j]==delta[j]*s);
                        n_s += 1;
                    }
                    ci = ci / 2;
                }
                
                bool is_s = true;
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
                if(!normal(n_s).is_zero()) {
                    is0 = true;
                    break;
                }
            }
            if(is0) break;
        }
        if(!is0) return lst{fun, exp, deltas};
        else return lst{};
    }, "SL", 0, true);
    
    FunExp.clear();
    FunExp.shrink_to_fit();
    for(auto item : sl_res) {
        lst fed = ex_to<lst>(item);
        if(fed.nops()<1) continue;
        FunExp.push_back(fed);
    }
    
    if(verb) cout << FunExp.size() << endl;
    if(FunExp.size()<1) IsZero = true;
}

void SD::RemoveDeltas() {
    if(IsZero) return;
    
    while(true) {
        bool exit = true;
        vector<lst> funexp = FunExp;
        FunExp.clear();
        FunExp.shrink_to_fit();
        for(auto fe : funexp) {
            if(fe.nops()<3) {
                FunExp.push_back(fe);
                continue;
            }
            auto xs = ex_to<lst>(fe.op(2).op(0));
            lst re_deltas;
            for(int i=1; i<fe.op(2).nops(); i++) {
                re_deltas.append(fe.op(2).op(i));
            }
            if(re_deltas.nops()>0) exit = false;
            
            if(xs.nops()<1) {
                FunExp.push_back(lst{fe.op(0), fe.op(1), re_deltas});
            } else {
                for(int i=0; i<xs.nops(); i++) {
                    auto xj = xs.op(i);
                    auto jInv = lstHelper::sum(xs).subs(xj==1);
                    
                    exmap repl;
                    for(int j=0; j<xs.nops(); j++) {
                        auto xxj = xs.op(j);
                        if(xxj != xj) repl[xxj] = xj*xxj;
                    }
                    
                    lst funs;
                    lst exps = ex_to<lst>(fe.op(1));
                    ex expns = 0;
                    for(int j=0; j<fe.op(0).nops(); j++) {
                        auto fun = fe.op(0).op(j);
                        fun = fun.subs(repl).normal();
                        if(!fun.is_polynomial(xj)) {
                            cerr << "xj: " << xj << endl;
                            cerr << "fun: " << fun << endl;
                            cerr << "fun is NOT polynormial of xj." << endl;
                            assert(false);
                        }
                        auto expn = expand(fun).degree(xj);
                        fun = pow(xj, -expn) * fun;
                        fun = normal(fun.subs(xj==1/xj));
                        fun = fun.subs(xj==jInv);
                        funs.append(fun);
                        expns += expn * exps.op(j);
                    }
                    
                    funs.append(jInv);
                    exps.append(ex(0)-xs.nops()-expns);
                    if(re_deltas.nops()>0) FunExp.push_back(lst{funs, exps, re_deltas});
                    else FunExp.push_back(lst{funs, exps});
                }
            }
        }
        if(exit) break;
    }
    
    XReOrders();
    Normalizes();
}

void SD::XEnd() {
    if(Verbose > 0) cout << now() << " - BiSection: " << FunExp.size() << " :> " << flush;
    vector<ex> funexps =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, FunExp, [&](auto &fe, auto rid) {
        lst para_res_lst;
        if(xSign(fe.op(0).op(1))==0) {
            auto fes = AutoEnd(fe);
            for(auto fei : fes) para_res_lst.append(fei);
        } else {
            para_res_lst.append(fe);
        }
        return para_res_lst;
    }, "f1", 0, !debug);

    FunExp.clear();
    FunExp.shrink_to_fit();
    for(auto &item : funexps) {
        for(auto &it : ex_to<lst>(item)) FunExp.push_back(ex_to<lst>(it));
    }
    if(Verbose > 0) cout << FunExp.size() << endl;
}

// after SDPrepares, Integrands can be expanded in ep safely.
void SD::SDPrepares() {
    if(IsZero) return;
    if(FunExp.size()<1) {
        IsZero = true;
        return;
    }
    if(SecDec==NULL) SecDec = new SecDecG();
    
    lst isyms = { ep, eps, vs, vz, iEpsilon };
    for(auto is : isyms) ParallelSymbols.append(is);
    ParallelSymbols.sort();
    ParallelSymbols.unique();
    
    MB();
    RemoveDeltas();
    if(CheckEnd) XEnd();
    
    Integrands.clear();
    Integrands.shrink_to_fit();
    
    auto fes = FunExp;
    FunExp.clear();
    FunExp.shrink_to_fit();
    for(auto &fe : fes) {
        bool to_add = true;
        for(auto item : fe.op(0)) {
            if(item.is_zero()) {
                to_add = false;
                break;
            }
        }
        if(to_add) FunExp.push_back(fe);
    }
    if(FunExp.size()<1) {
        IsZero = true;
        return;
    }
    
    if(Verbose > 0) cout << now() << " - SDPrepares: ..." << endl << flush;
    
    vector<ex> sd_res =
    GiNaC_Parallel(ParallelProcess, ParallelSymbols, FunExp, [&](auto &fe, auto rid) {
        // return a lst, element pattern: { {{x1,n1}, {x2,n2}, ...}, {{e1, n1},{e2,n2}, ...} }.
        lst para_res_lst;
        auto xns_pns = DS(fe);
        for(auto const &item : xns_pns) {

            // take z-poles
            if(item.has(vz)) {
                auto ct = item.op(1).op(0).op(0);
                ct = ct.subs(lst{ CT(wild())==wild(),FTX(wild(1),wild(2))==1 }).subs(pow(vs,vz)==1);
                int sNN = sN - vsRank(ct.subs(pow(vs,vz)==1));
            
                lst zpols;
                // poles from Gamma(-z)
                for(int vn=0; vn<=sNN; vn++) {
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
                            if(zpn>sNN) break;
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
                } else {
                    assert(false);
                }
            }
            if(tmp!=1) pns.append(lst{tmp, 1});
            ibp_in_vec.push_back(lst{xns, pns});
        }
    }
    if(Verbose > 1) cout << "  \\--" << WHITE << "Maximum x^n: " << ex(0)-min_expn << " + 1" << RESET << endl << flush;

    int pn = 0;
    vector<ex> ibp_res_vec;
    while(ibp_in_vec.size()>0) {
        pn++;
        ostringstream spn;
        spn << "IBP-" << pn;
        vector<ex> ibp_res =
        GiNaC_Parallel(ParallelProcess, ParallelSymbols, ibp_in_vec, [&](auto &xns_pns, auto rid) {
            // return lst
            // {0, element} for input with pole reached and doing nothing
            // {1, {element, ...}} for input whth pole NOT reached
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
                    
                    for(int i=0; i<pns.nops(); i++) {
                        lst xns3 = ex_to<lst>(xns);
                        xns3.let_op(n).let_op(1) = xn.op(1)+1;
                        
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
                    return lst{1, xns_pns_lst};
                }
            }
            return lst{0, xns_pns };

        }, spn.str().c_str(), Verbose, true);
    
        ibp_in_vec.clear();
        ibp_in_vec.shrink_to_fit();
        for(auto &ii : ibp_res) {
            auto check = ii.op(0);
            if(check>0) {
                auto items = ii.op(1);
                for(auto &it : ex_to<lst>(items)) ibp_in_vec.push_back(it);
            } else {
                auto item = ii.op(1);
                ex expr = 1;
                for(auto pn : item.op(1)) expr *= pow(pn.op(0), pn.op(1));
                ibp_res_vec.push_back(lst{ item.op(0), expr });
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
                auto ct2 = vc * ct;
                int ctN = epRank(ct2);
                tmp = mma_series(tmp, ep, epN-ctN);
                for(int di=tmp.ldegree(ep); (di<=tmp.degree(ep) && di<=epN-ctN); di++) {
                    auto intg = tmp.coeff(ep, di);
                    assert(!intg.has(ep));
                    auto pref = mma_series(ct2, ep, epN-di);
                    if(pref.has(vs)) pref = mma_series(pref, vs, sN);
                    //if(use_CCF) intg = collect_common_factors(intg);
                    para_res_lst.append(lst{pref * pow(ep, di), intg});
                }
            } else {
                auto sct = vc * ct;
                int sctN = epsRank(sct);
                ex stmp = mma_series(tmp, eps, epsN-sctN);
                for(int sdi=stmp.ldegree(eps); (sdi<=stmp.degree(eps) && sdi<=epsN-sctN); sdi++) {
                    tmp = stmp.coeff(eps, sdi);
                    //if(use_CCF) tmp = collect_common_factors(tmp);
                    assert(!tmp.has(eps));
                    auto ct2 = mma_series(sct, eps, epsN-sdi);
                    int ctN = epRank(ct2);
                    tmp = mma_series(tmp, ep, epN-ctN);
                    for(int di=tmp.ldegree(ep); (di<=tmp.degree(ep) && di<=epN-ctN); di++) {
                        auto intg = tmp.coeff(ep, di);
                        assert(!intg.has(ep));
                        auto pref = mma_series(ct2, ep, epN-di);
                        if(pref.has(vs)) pref = mma_series(pref, vs, sN);
                        //if(use_CCF) intg = collect_common_factors(intg);
                        para_res_lst.append(lst{pref * pow(eps, sdi) * pow(ep, di), intg});
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

}
