/**
 * @file
 * @brief Basic Functions for SecDec
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */
 
#include "SD.h"
#include <math.h>
#include <cmath>
#include<algorithm>

namespace HepLib::SD {

    bool SecDecBase::VerifySD(vector<exmap> map_vec, bool quick) {
        lst xs;
        exmap nxs;
        for(auto kv : map_vec[0]) {
            auto xi = kv.first;
            if(is_zero(xi-x(-1))) continue;
            xs.append(xi);
            if(quick) nxs[xi] = xi.subs(x(w)==ex(w+1)/ex(w+13));
            else nxs[xi];
        }
        ex chk_total = 1;
        for(auto xi : xs) chk_total *= 1/(nxs[xi]+1);
        
        ex int_total = 0;
        for(auto imap : map_vec) {
            ex intg = 1;
            for(auto kv : imap) {
                auto xi = kv.first;
                auto nxi = kv.first.subs(x(w)==(13+w));
                if(is_zero(xi-x(-1))) intg *= kv.second;
                else intg *= pow(kv.second, nxs[kv.first]);
            }
            while(true) {
                auto tmp = intg;
                tmp = tmp.subs(w0*pow(w1*w2,w3)==w0*pow(w1,w3)*pow(w2,w3));
                tmp = xyz_pow_simplify(tmp);
                if(is_zero(tmp-intg)) break;
                intg = tmp;
            }
            intg = intg.subs(pow(y(w1),w2)==1/(w2+1)).subs(y(w)==1/ex(2));
            int_total += intg;
        }
        int_total = int_total.normal();
        return normal(chk_total-int_total).is_zero();
    }
    bool SecDec::VerifySD(vector<exmap> map_vec, bool quick) {
        return SecDecBase::VerifySD(map_vec, quick);
    }

    vector<exmap> SecDecBase::x2y(const lst &in_xpols) {
        if(in_xpols.has(y(w))) throw Error("SecDecBase::x2y: y(w) found @ " + ex2str(in_xpols));
        lst xpols_lst = in_xpols;
        while(true) {
            ex xpols = 1;
            for(auto item : xpols_lst) {
                exmap pmap;
                item = item.to_rational(pmap);
                xpols *= factor(item).subs(pmap);
            }
            if(!is_a<mul>(xpols)) xpols = lst{ xpols };
            
            xpols_lst.remove_all();
            for(auto item : xpols) {
                if(item.match(x(w))) continue;
                else if(item.match(pow(x(w1),w2))) continue;
                else if(item.match(pow(w1,w2))) {
                    xpols_lst.append(item.op(0));
                } else {
                    xpols_lst.append(item);
                }
            }
            xpols_lst.sort();
            xpols_lst.unique();
            
            ex xpols2 = 1;
            for(auto item : xpols_lst) xpols2 *= item;
            if(is_a<lst>(xpols)) xpols = xpols.op(0);
            if(is_zero(xpols2-xpols)) return x2y(xpols);
        }
        throw Error("SecDecBase::x2y: unexpected region reached.");
        return vector<exmap>();
    }

    /*-----------------------------------------------------*/
    /*                   CheckAtEnd                        */
    /*-----------------------------------------------------*/
    bool SecDec::IsBad(ex f, vector<exmap> vmap) {
        for(auto &vi : vmap) {
            auto ft = f.subs(vi);
            auto xs_tmp = get_x_from(ft);
            auto ys_tmp = get_y_from(ft);
            int ysn = ys_tmp.size();
            vector<int> y_free_indexs;
            for(int j=0; j<ys_tmp.size()+10; j++) {
                if(!ft.has(y(j))) y_free_indexs.push_back(j);
            }
            if(xs_tmp.size()>y_free_indexs.size()) throw Error("IsBad: xsize>ysize");
            for(int i=0; i<xs_tmp.size(); i++) {
                vi[xs_tmp[i]] = y(y_free_indexs[i]);
            }
            if(xs_tmp.size()>0) ft = ft.subs(vi);

            // need collect_common_factors
            ft = collect_common_factors(mma_expand(ft,y(w)));
            if(is_exactly_a<mul>(ft)) {
                ex ret = 1;
                for (auto item : ft) {
                    if( !(item.match(y(w)) || item.match(pow(y(w), w1))) ) {
                        ret *= item;
                    }
                }
                ft = ret;
            } else if ( ft.match(y(w)) || ft.match(pow(y(w), w1)) ) {
                ft = 1;
            } // else ft not changed

            ys_tmp = get_y_from(ft);
            for(int pi=1; pi<std::pow(2, ys_tmp.size()); pi++) {
                lst yrepl;
                int ci = pi;
                for(int i=0; i<ys_tmp.size(); i++) {
                    yrepl.append(ys_tmp[i] == ((ci%2)==1 ? 1 : 0));
                    ci /= 2;
                }
                if(normal(ft.subs(yrepl)).is_zero()) return true;
            }
        }
        return false;
    }

    vector<ex> SecDec::AutoEnd(ex po_ex) {
        if(po_ex.nops()>2) throw Error("AutoEnd: Deltas found @ " + ex2str(po_ex));
        lst const exlist = ex_to<lst>(po_ex.op(1));
        if(!(exlist.op(0)-1).is_zero()) throw Error("AutoEnd: (!(exlist.op(0)-1).is_zero())");
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
            for(auto polist : polists) {
                lst sdList;
                for(int i=0; i<polist.nops(); i++) {
                    auto tmp = polist.op(i);
                    auto ntmp = NN(exlist.op(i));
                    if(!tmp.subs(lst{x(w)==0, y(w)==0}).normal().is_zero()) continue;
                    if( (!tmp.has(x(w)) && !tmp.has(y(w))) || (is_a<numeric>(ntmp) && ntmp>0) ) continue;
                    sdList.append(tmp);
                }
                vector<exmap> vmap = SecDec->x2y(sdList);
                
                for(int ni=0; ni<polist.nops(); ni++) {
                    auto po = polist.op(ni);
                    auto expo = exlist.op(ni); // note that expo may decrease when using diff
                    if(expo.info(info_flags::posint)) continue;
                    if(IsBad(po, vmap)) {
                        OK = false;
                        break;
                    }
                }
                if(!OK) break;
            }

            if(OK) {
                vector<ex> res;
                for(auto item : polists) res.push_back(lst{ex_to<lst>(item), exlist});
                return res;
            }
        }}
        
        cerr << ErrColor << "polynormial list: " << po_ex.op(0) << RESET << endl;
        throw Error("AutoEnd Failed @ ALL possible bisections!");
        return vector<ex>();
    }

    /*-----------------------------------------------------*/
    /*                       SD                            */
    /*-----------------------------------------------------*/
    //return a lst, element pattern: { {{x1,n1}, {x2,n2}, ...}, {{e1, n1},{e2,n2}, ...} }
    //e1 is a const term, e2 still the F-term
    vector<ex> SecDec::DS(const ex po_ex) {
        // 1st element in input polist is the constant term, guess NOT necessary
        // 2nd element in input polist is the F-term, required!
        lst const polist = ex_to<lst>(po_ex.op(0));
        lst const exlist = ex_to<lst>(po_ex.op(1));
        lst sdList;
        for(int i=0; i<polist.nops(); i++) {
            auto tmp = polist.op(i);
            auto ntmp = exlist.op(i);
            if(!tmp.subs(lst{x(w)==0, y(w)==0}).normal().is_zero()) continue;
            ntmp = NN(ntmp);
            if( (!tmp.has(x(w)) && !tmp.has(y(w))) || (is_a<numeric>(ntmp) && ntmp>0) ) continue;
            sdList.append(tmp);
        }

        vector<exmap> vmap = SecDec->x2y(sdList);
        if(!VerifySD(vmap)) {
            for(auto vi : vmap) cout << vi << endl;
            throw Error("VerifySD Failed!");
        }

        vector<ex> sd;
        for(auto vi : vmap) {
            auto ypolist = polist.subs(vi);
            auto xs_tmp = get_x_from(ypolist);
            auto ys_tmp = get_y_from(ypolist);
            int ysn = ys_tmp.size();
            vector<int> y_free_indexs;
            for(int j=0; j<xs_tmp.size()+ys_tmp.size()+10; j++) {
                if(!ypolist.has(y(j))) y_free_indexs.push_back(j);
            }
            if(xs_tmp.size()>y_free_indexs.size()) throw Error("DS: xsize>ysize");
            for(int i=0; i<xs_tmp.size(); i++) {
                vi[xs_tmp[i]] = y(y_free_indexs[i]);
            }
            if(xs_tmp.size()>0) ypolist = ypolist.subs(vi);
            if(ypolist.has(x(w))) throw Error("DS: x(w) found @ " + ex2str(ypolist));

            // need collect_common_factors
            auto ft = collect_common_factors(mma_expand(ypolist.op(1),y(w)));
            
            ex ct = 1, fsgin = 1;
            if(is_a<mul>(ft)) {
                ex ret = 1;
                for (auto item : ft) {
                    if( !(item.match(y(w)) || item.match(pow(y(w), w1))) ) {
                        ret *= item;
                    }
                }
                ft = ret;
            } else if( ft.match(y(w)) || ft.match(pow(y(w), w1)) ) {
                ft = 1;
            } // else ft not changed
            bool hasWRA = false;
            if(ft.has(WRA(w))) {
                ft = 1;
                hasWRA = true;
            }

            bool need_contour_deformation = ft.has(PL(w));
            if(ft.has(y(w)) && !need_contour_deformation) {
                auto tmp = ft.subs(nReplacements).subs(lst{
                    CV(w1,w2)==w2, ep==ex(1)/111, eps==ex(1)/1111
                }).expand();
                if(is_a<add>(tmp)) {
                    need_contour_deformation = false;
                    auto first = tmp.op(0).subs(y(w)==1);
                    for(auto item : tmp) {
                        auto chk = NN(item.subs(y(w)==1)*first);
                        if(!is_a<numeric>(chk)) {
                            throw Error("DS: Not a number: " + ex2str(item.subs(y(w)==1)*first));
                        }
                        if(chk<0) {
                            need_contour_deformation = true;
                            break;
                        }
                    }
                }
            }

            if(!need_contour_deformation) {
                auto tmp = NN(ft.subs(y(w)==1).subs(nReplacements).subs(lst{
                    CV(w1,w2)==w2, ep==ex(1)/111, eps==ex(1)/1111
                }));
                
                if(!is_a<numeric>(tmp)) throw Error("DS: NOT a numeric with " + ex2str(tmp));
                if(tmp<0) {
                    ct = exp(-I * Pi * exlist.op(1));
                    fsgin = -1;
                }
                ft = 1;
            }

            exmap ymol;
            
            // need collect_common_factors
            auto det = collect_common_factors(vi[x(-1)]);
            if(is_a<add>(det)) throw Error("DS: det is add " + ex2str(det));
            auto ys = get_xy_from(det);
            ex det1 = 1;
            for(int i=0; i<ys.size(); i++) {
                auto ndeg = det.degree(ys[i]);
                if(ndeg!=det.ldegree(ys[i])) throw Error("DS: det is not a monomial with det = "+ex2str(det));
                ymol[ys[i]] = ymol[ys[i]] + ndeg;
                det1 *= pow(ys[i], ndeg);
            }
            if(!(is_a<numeric>(det/det1) && ex_to<numeric>(det/det1).is_integer())) {
                throw Error("DS: det=" + ex2str(det) + ", det1=" + ex2str(det1));
            }
            
            lst ftxlst;
            for(auto xi : get_xy_from(ft)) ftxlst.append(xi);
            
            lst pol_exp_lst;
            pol_exp_lst.append(lst{ subs(FTX(ft,ftxlst)*CT(ct*det/det1), y(w)==x(w)), 1 });

            for(int i=0; i<ypolist.nops(); i++) {
                auto tmp = ypolist.op(i);
                auto nexp = exlist.op(i).normal();
                bool nchk = (!is_a<numeric>(nexp) || (is_a<numeric>(nexp) && nexp<0));
                if(nexp.has(x(w)) || nexp.has(y(w))) throw Error("DS: x or y found in exp: "+ex2str(nexp));
                
                if(tmp.has(y(w))) {
                    lst ys;
                    for(auto yi : get_y_from(tmp)) ys.append(yi);
                    if(!tmp.is_polynomial(ys)) throw Error("DS: NOT a polynomial with "+ex2str(tmp));
                }
        
                if(i==1 && fsgin<0) tmp = ex(0)-tmp; // check complete negtive F-term
                auto tmp1 = tmp;
                while(true) {
                    tmp=tmp1.subs(pow(y(w1)*w2, w3)==pow(y(w1),w3) * pow(w2, w3));
                    tmp = tmp.subs(pow(pow(y(w1),w2),w3)==pow(y(w1),w2*w3));
                    tmp = tmp.subs(pow(sqrt(y(w1)),w2)==pow(y(w1),w2/2));
                    tmp = tmp.subs(sqrt(pow(y(w1),w2))==pow(y(w1),w2/2));
                    if((tmp1-tmp).is_zero()) break;
                    tmp1 = tmp;
                }
                
                // need collect_common_factors
                if(tmp.has(y(w))) tmp = collect_common_factors(mma_expand(tmp,y(w)));

                lst tmps;
                if(is_exactly_a<mul>(tmp)) {
                    for (auto item : tmp) tmps.append(item);
                } else {
                    tmps.append(tmp);
                }
                
                ex rem = 1;
                ex ct = 1;
                for (auto item : tmps) {
                    if(item.match(y(w))) {
                        auto yi = item;
                        ymol[yi] = ymol[yi] + nexp;
                    } else if(item.match(pow(y(w1), w2))) {
                        auto yi = item.op(0);
                        ymol[yi] = ymol[yi] + item.op(1) * nexp;
                    } else if(!item.has(y(w)) && !item.has(x(w))) {
                        if(is_a<numeric>(nexp) && ex_to<numeric>(nexp).is_integer()) {
                            ct *= item;
                        } else if(!item.has(PL(w)) && !item.has(WRA(w))) {
                            auto tr = NN(item.subs(nReplacements).subs(lst{
                                CV(w1,w2)==w2, ep==ex(1)/111, eps==ex(1)/1111
                            }));
                            if(!is_a<numeric>(tr)) {
                                throw Error("DS: not numeric - item: " + ex2str(tr) + " ; " + ex2str(item));
                            }
                            auto nitem = ex_to<numeric>(tr);
                            if( nitem.is_real() && nitem<0 ) {
                                ct *= ex(-1)*item;
                                rem *= ex(-1);
                            } else {
                                ct *= item;
                            }
                        } else {
                            rem *= item;
                        }
                    } else {
                        if(nchk && item.subs(lst{y(w)==0,iEpsilon==0}).normal().is_zero()) {
                            throw Error("DS: zero item: " + ex2str(item)  + " and exlist.op(i) = " + ex2str(nexp));
                        }
                        rem *= item;
                    }
                }

                ex pnp = rem-(i==1 && (ft!=1 || hasWRA) ? iEpsilon : ex(0));
                pnp = pnp.subs(y(w)==x(w));
                
                pol_exp_lst.append(lst{pnp, nexp});
                pol_exp_lst.let_op(0) = pol_exp_lst.op(0).subs(CT(w)==CT(w*pow(ct,nexp))).subs(CT(0)==0);
            }

            lst x_n_lst;
            for(auto & kv : ymol) {
                auto k = kv.first.subs(y(w)==x(w));
                auto v = kv.second;
                if(is_a<numeric>(v)) {
                    auto nv = ex_to<numeric>(v);
                    if(nv<=-1) {
                        throw Error("DS: " + ex2str(k) + "^(" + ex2str(nv) + ") found, more regularization needed!");
                    }
                }
                x_n_lst.append(lst{k, v});
            }
            x_n_lst.sort();
            sd.push_back(lst{x_n_lst, pol_exp_lst});
        }
        
        return sd;
    }

    // 1st element in output is the constant term
    // 2nd element in both input and output is the F-term
    lst SecDec::Normalize(const ex &input) {
        ex const_term = 1;
        lst plst, nlst;
        lst in_plst = ex_to<lst>(input.op(0));
        lst in_nlst = ex_to<lst>(input.op(1));
        for(int i=0; i<in_plst.nops(); i++) {
            if(i!=1 && (in_nlst.op(i).is_zero() || in_plst.op(i)==ex(1))) continue;
            if(i!=1 && !in_plst.op(i).has(x(w)) && !in_plst.op(i).has(y(w)) && !in_plst.op(i).has(z(w))) {
                if(in_nlst.op(i).info(info_flags::integer)) const_term *= pow(in_plst.op(i), in_nlst.op(i));
                else const_term *= exp(log(in_plst.op(i)) * in_nlst.op(i));
            } else {
                auto ptmp = collect_common_factors(mma_expand(in_plst.op(i),{x(w),y(w),z(w)}));
                auto ntmp = in_nlst.op(i);
                if(is_a<mul>(ptmp)) {
                    ex tmul = 1;
                    for(int j=0; j<ptmp.nops(); j++) {
                        auto tmp = ptmp.op(j);
                        if(!tmp.has(x(w)) && !tmp.has(y(w)) && !tmp.has(z(w))) { // constant terms
                            if(ntmp.info(info_flags::integer)) {
                                const_term *=  pow(tmp,ntmp);
                            } else if((tmp-vs).is_zero() || tmp.match(pow(vs,w))) {
                                const_term *=  pow(tmp,ntmp);
                            } else if(!tmp.has(PL(w)) && !tmp.has(vs) && !tmp.has(WRA(w))) {
                                auto tr = NN(tmp.subs(nReplacements).subs(lst{
                                    CV(w1,w2)==w2, ep==ex(1)/111, eps==ex(1)/1111
                                }));
                                if(!is_a<numeric>(tr) || !tr.info(info_flags::real)) {
                                    cerr << "tmp: " << tmp << endl;
                                    cerr << "tr: " << tr << endl;
                                    cerr << "nReplacements: " << nReplacements << endl;
                                    throw Error("Normalize: tr is NOT numeric with nReplacements.");
                                }
                                
                                if(tr>0) {
                                    const_term *=  pow(tmp,ntmp);
                                } else {
                                    const_term *= pow(-tmp,ntmp);
                                    tmul = ex(0) - tmul;
                                }
                            } else {
                                tmul *= tmp;
                            }
                        } else if(tmp.match(pow(w1,w2)) && xPositive(tmp.op(0))) {
                            plst.append(tmp.op(0));
                            nlst.append(ntmp * tmp.op(1));
                        } else if(xSign(tmp)!=0 || xSign(tmp.subs(lst{y(w)==x(w),z(w)==x(w)}))!=0) {
                            if(tmp.subs(lst{x(w)==1, y(w)==1, z(w)==1})>0) {
                                plst.append(tmp);
                                nlst.append(ntmp);
                            } else {
                                plst.append(ex(0)-tmp);
                                nlst.append(ntmp);
                                tmul = ex(0) - tmul;
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
                if(inp[plst[i]]==1) continue;
                plst_comb.append(plst[i]);
                nlst_comb.append(np[plst[i]]);
                inp[plst[i]] = 1;
            }
        }
        
        if(nlst_comb.op(0)!=1) {
            if(nlst_comb.op(0).info(info_flags::integer)) plst_comb.let_op(0) = pow(plst_comb.op(0),nlst_comb.op(0));
            else plst_comb.let_op(0) = exp(log(plst_comb.op(0))*nlst_comb.op(0));
            nlst_comb.let_op(0) = 1;
        }
        
        return (input.nops()>2) ? lst{plst_comb, nlst_comb, input.op(2)} : lst{plst_comb, nlst_comb};
    }

    /*-----------------------------------------------------*/
    /*               's Funtions in SD                     */
    /*-----------------------------------------------------*/
    // working with or without Deltas
    void SecDec::XReOrders() {
        if(IsZero || !use_XReOrders) return;
        if(Integrands.size()<1) {
            for(auto &fe : FunExp) {
                auto xs = get_x_from(fe);
                lst x2y;
                for(int i=0; i<xs.size(); i++) x2y.append(xs[i]==y(i));
                fe = ex_to<lst>(subs(fe, x2y).subs(y(w)==x(w)));
            }
        } else {
            for(auto &vint : Integrands) {
                auto xs = get_x_from(vint);
                lst x2y;
                for(int i=0; i<xs.size(); i++) x2y.append(xs[i]==y(i));
                vint = vint.subs(x2y).subs(y(w)==x(w));
            }
        }
    }

    // working with or without Deltas
    void SecDec::Normalizes() {
        for(int ri=0; ri<2; ri++) { // run twice, needs to check in more details
            if(IsZero) return;
            
            vector<ex> funexp;
            for(auto fe : FunExp) {
                funexp.push_back(Normalize(fe));
            }
            FunExp.clear();
            FunExp.shrink_to_fit();
            
            exmap fn;
            for(auto fe : funexp) {
                ex key = 1;
                if(fe.nops()>2) key = iWF(fe.op(2)); // deltas
                if(fe.op(1).op(0)!=1) {
                    cout << fe << endl;
                    throw Error("Normalizes: fe.op(0).op(0) is NOT 1.");
                }
                for(int i=1; i<fe.op(0).nops(); i++) key *= pow(fe.op(0).op(i), fe.op(1).op(i));
                fn[key] += fe.op(0).op(0);
            }
            
            map<ex,int,ex_is_less> ifn;
            for(auto fe : funexp) {
                ex key = 1;
                if(fe.nops()>2) key = iWF(fe.op(2));
                for(int i=1; i<fe.op(0).nops(); i++) key *= pow(fe.op(0).op(i), fe.op(1).op(i));
                if(ifn[key]==1) continue;
                lst funs, exps;
                funs.append(fn[key]);
                exps.append(1);
                for(int i=1; i<fe.op(0).nops(); i++) {
                    funs.append(fe.op(0).op(i));
                    exps.append(fe.op(1).op(i));
                }
                if(fe.nops()>2) FunExp.push_back(lst{funs, exps, fe.op(2)});
                else FunExp.push_back(lst{funs, exps});
                ifn[key] = 1;
            }
        }
    }

    // working with or without Deltas
    void SecDec::XTogethers() {
        vector<ex> funexp;
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
                if(expi.info(info_flags::nonnegint)) {
                    rem *= pow(fe.op(0).op(i), expi);
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
    void SecDec::XExpands() {
        vector<ex> funexp;
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
                if(expi.info(info_flags::nonnegint)) {
                    rem *= pow(fe.op(0).op(i), expi);
                } else {
                    fun.append(fe.op(0).op(i));
                    exp.append(fe.op(1).op(i));
                }
            }
            
            rem = mma_collect(rem, x(w));
            
            lst rem_lst;
            if(is_a<add>(rem)) {
                for(auto item : rem) rem_lst.append(item);
            } else {
                rem_lst.append(rem);
            }
            
            for(auto item : rem_lst) {
                auto fe2 = fe;
                let_op_append(fe2, 0, item);
                let_op_append(fe2, 1, 1);
                fe2.let_op(0) = fun;
                fe2.let_op(1) = exp;
                FunExp.push_back(fe2);
            }
        }
    }

    // Section 2.1 @ https://arxiv.org/pdf/1712.04441.pdf
    // also refers to Feng/Thinking.pdf
    void SecDec::Scalelesses() {
        if(IsZero) return;
        if(Verbose > 2) cout << "  \\--Scaleless: " << FunExp.size() << " :> " << flush;
        
        auto verb = Verbose; Verbose = 0;
        auto sl_res =
        GiNaC_Parallel(FunExp.size(), [&](int idx)->ex {
            auto funexp = FunExp[idx];
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
                        if(exp.op(j).info(info_flags::nonnegint)) continue;
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
        }, "SL", true);
        Verbose = verb;
        
        FunExp.clear();
        FunExp.shrink_to_fit();
        for(auto item : sl_res) {
            lst fed = ex_to<lst>(item);
            if(fed.nops()<1) continue;
            FunExp.push_back(fed);
        }
        
        if(Verbose > 2) cout << FunExp.size() << endl;
        if(FunExp.size()<1) IsZero = true;
    }

    void SecDec::RemoveDeltas() {
        if(IsZero) return;
        
        while(true) {
            bool exit = true;
            vector<ex> funexp = FunExp;
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
                        ex jInv = 0;
                        for(auto xi : xs) jInv += xi;
                        jInv = jInv.subs(xj==1);
                        
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
                                cerr << funexp << endl;
                                throw Error("RemoveDeltas: fun is NOT polynormial of xj.");
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

    void SecDec::XEnd() {
        if(Verbose > 2) cout << "  \\--BiSection: " << FunExp.size() << " :> " << flush;
        auto verb = Verbose; Verbose=0;
        auto funexps =
        GiNaC_Parallel(FunExp.size(), [&](int idx)->ex {
            auto fe = FunExp[idx];
            lst para_res_lst;
            if(xSign(fe.op(0).op(1))==0 && !fe.has(WRA(w))) {
                auto fes = AutoEnd(fe);
                for(auto fei : fes) para_res_lst.append(fei);
            } else {
                para_res_lst.append(fe);
            }
            return para_res_lst;
        }, "BiSec", !debug);
        Verbose = verb;

        FunExp.clear();
        FunExp.shrink_to_fit();
        for(auto &item : funexps) {
            for(auto &it : ex_to<lst>(item)) FunExp.push_back(ex_to<lst>(it));
        }
        if(Verbose > 2) cout << FunExp.size() << endl;
    }

    // after SDPrepares, Integrands can be expanded in ep safely.
    void SecDec::SDPrepares() {
        if(IsZero) return;
        if(FunExp.size()<1) {
            IsZero = true;
            return;
        }
        if(SecDec==NULL) SecDec = new SecDecG();
        
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
        
        if(Verbose > 1) cout << Color_HighLight << "  SDPrepares @ " << now() << RESET << endl;
        auto sd_res =
        GiNaC_Parallel(FunExp.size(), [&](int idx)->ex {
            // return a lst, element pattern: { {{x1,n1}, {x2,n2}, ...}, {{e1, n1},{e2,n2}, ...} }.
            auto fe = FunExp[idx];
            auto xns_pns = DS(fe);
            lst para_res_lst;
            for(auto const &item : xns_pns) {
        
                // take z-poles
                if(item.has(vz)) {
                    auto ct = item.op(1).op(0).op(0);
                    ct = ct.subs(lst{ CT(w)==w,FTX(w1,w2)==1 }).subs(pow(vs,vz)==1);
                    int sNN = sN - vsRank(ct.subs(pow(vs,vz)==1));

                    lst zpols;
                    // poles from Gamma(-z)
                    for(int vn=0; vn<=sNN; vn++) zpols.append(vn);

                    // poles from xi^{c1*z+c0} with c1<0
                    for(auto xn : item.op(0)) {
                        auto xn_op1 = xn.op(1).normal().expand();
                        ex c1 = xn_op1.coeff(vz);
                        ex c0 = xn_op1.subs(vz==0);
                        
                        if(!is_a<numeric>(c1)) {
                            cerr << ErrColor << "SDPrepares: c1 is not a number: " << c1 << RESET << endl;
                            exit(1);
                        }
                        
                        if(c1<0) {
                            int pxn = -1;
                            while(true) {
                                ex zp = (pxn-c0)/c1;
                                ex zpn = zp.subs(lst{eps==0,ep==0,epz==0});
                                if(!is_a<numeric>(zpn)) {
                                    cerr << ErrColor << "SDPrepares: zpn is not a number: " << zpn << RESET << endl;
                                    exit(1);
                                }
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
                        auto item2 = item.subs(vz==ss+zp).subs(ss==vz);
                        para_res_lst.append(item2);
                    }
                } else {
                    para_res_lst.append(item);
                }
            }
            return para_res_lst;
        }, "SD", true);
        
        ex min_expn = 1, min_expn2 = 10;
        vector<ex> ibp_in_vec;
        for(auto &item : sd_res) {
            for(auto &it : ex_to<lst>(item)) {
                ex expn = 0;
                for(auto xn : it.op(0)) {
                    ex nxn = xn.op(1).subs(lst{ep==0, eps==0, vz==0});
                    if(nxn<-1) expn += nxn+1;
                    if(min_expn2>nxn) min_expn2 = nxn+1;
                }
                if(expn < min_expn) min_expn = expn;
                
                lst xns = ex_to<lst>(it.op(0));
                lst pns;
                for(int i=0; i<it.op(1).nops(); i++) {
                    lst pn = ex_to<lst>(it.op(1).op(i));
                    if(i<2) {
                        pns.append(pn);
                        continue;
                    }
                    if(pn.op(0).is_equal(1)) continue;
                    if(is_zero(pn.op(1))) {
                        if(is_zero(pn.op(0))) throw Error("SDPrepares: 0^0 found.");
                        continue;
                    }
                    pns.append(pn);
                }
                ibp_in_vec.push_back(lst{xns, pns});
            }
        }

        if(Verbose > 1) cout << "  \\--" << Color_HighLight << "Maximum x^-n: All(" << ex(0)-min_expn << "+1X), Max(" << (ex(0)-min_expn2) << "+1)" << RESET << endl;

        int pn = 0;
        vector<ex> ibp_res_vec;
        while(ibp_in_vec.size()>0) {
            pn++;
            ostringstream spn;
            spn << "IBP-" << (pn-1);
            auto ibp_res =
            GiNaC_Parallel(ibp_in_vec.size(), [&](int idx)->ex {
                // return lst
                // {0, element} for input with pole reached and doing nothing
                // {1, {element, ...}} for input whth pole NOT reached
                // element pattern still as { {{x1,n1}, {x2,n2}, ...}, {{e1, n1},{e2,n2}, ...} }
                
                auto xns_pns = ibp_in_vec[idx];
                
                auto xns = xns_pns.op(0);
                auto pns = xns_pns.op(1);
                
                exset fts;
                pns.op(0).find(FTX(w1,w2), fts);
                bool noFT = (fts.size()==1) && ( is_zero((*(fts.begin())).op(0)-1) );
                
                ex pole_requested = -1;
                if(noFT || PoleRequested > -1) pole_requested = PoleRequested;
                
                for(int n=0; n<xns.nops(); n++) {
                    ex xn = xns.op(n);
                    auto expn = xn.op(1).subs(lst{eps==0,ep==0,vz==0,epz==0}).normal();
                    if(!is_a<numeric>(expn)) throw Error("SDPrepares: expn NOT numeric: " + ex2str(expn));

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
                            
                            ex tmp = ex(0)-pns.op(i).op(1)*mma_diff(pns.op(i).op(0),xx);
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
                            if(tz) tmp = collect_common_factors(mma_expand(tmp,x(w)));
                            if(tmp.is_zero()) continue;
                            
                            if(tz && is_a<mul>(tmp)) {
                                ex rem = 1;
                                for(auto ii : tmp) {
                                    if(ii.match(pow(x(w), w2))) {
                                        bool t = true;
                                        for(int ij=0; ij<xns3.nops(); ij++) {
                                            if(xns3.op(ij).op(0)==ii.op(0)) {
                                                xns3.let_op(ij).let_op(1) += ii.op(1);
                                                t = false;
                                                break;
                                            }
                                        }
                                        if(t) xns3.append(lst{ii.op(0), ii.op(1)});
                                    } else if(ii.match(x(w))) {
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
                            if(is_zero(pns.op(i).op(1)-1)) {
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

            }, spn.str().c_str(), true);
        
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
                    for(auto pn : item.op(1)) {
                        if(pn.op(0).has(CT(w)) || pn.op(0).has(FTX(w1,w2))) {
                            if(pn.op(1)!=1) throw Error("SDPrepares: exponent of CT is NOT 1.");
                            expr *= pn.op(0);
                        } else expr *= exp(log(pn.op(0)) * pn.op(1));
                    }
                    expr = exp_simplify(expr);
                    ibp_res_vec.push_back(lst{ item.op(0), expr });
                }
            }
        }

        auto res =
        GiNaC_Parallel(ibp_res_vec.size(), [&](int idx)->ex {

            // return single element in which ep/eps can be expanded safely.
            auto xns_expr = ibp_res_vec[idx];
            lst para_res_lst;
            auto xns = xns_expr.op(0);
            auto expr = xns_expr.op(1);
            lst exprs = { expr };
            symbol dx;
            for(auto xn : xns) {
                auto expn = xn.op(1).subs(lst{eps==0,ep==0,vz==0,epz==0}).normal();
                if(!is_a<numeric>(expn)) throw Error("SDPrepares: Not a number with expn = " + ex2str(expn));
                
                lst exprs2;
                for(auto it : exprs) {
                    ex rem = pow(xn.op(0), xn.op(1)) * it;
                    if(ex_to<numeric>(expn)<=-1) {
                        ex dit = it;
                        ex dit0 = dit.subs(xn.op(0)==0);
                        ex ifact = 1;
                        if(!is_zero(dit0)) {
                            rem -= pow(xn.op(0), xn.op(1)) * dit0 / ifact;
                            exprs2.append(dit0/(xn.op(1)+1)/ifact);
                        }
                        for(int i=1; i+expn<0; i++) {
                            dit = mma_diff(dit, xn.op(0));
                            if(is_zero(dit)) break;
                            dit0 = dit.subs(xn.op(0)==0);
                            ifact *= i;
                            if(!is_zero(dit0)) {
                                rem -= pow(xn.op(0), xn.op(1)+i) * dit0 / ifact;
                                exprs2.append(dit0/(xn.op(1)+i+1)/ifact);
                            }
                        }
                    }
                    exprs2.append(rem);
                }
                exprs = exprs2;
            }

            for(auto const &it : exprs) {
                if(it.is_zero()) continue;
                auto xs = get_x_from(it);
                lst x2y;
                for(int i=0; i<xs.size(); i++) x2y.append(xs[i]==y(i));
                para_res_lst.append(it.subs(x2y).subs(y(w)==x(w)));
            }

            return para_res_lst;
        }, "Taylor", true);
        
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
            GiNaC_Parallel(ints.size(), [&](int idx)->ex {
                auto item = ints[idx];
                ex it = item;
                if(it.has(vz)) {
                    exset cts;
                    it.find(CT(w),cts);
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
            }, "zResidue", true);
        } else {
            Integrands = ints;
        }

    }

    void SecDec::EpsEpExpands() {
        if(IsZero) return;
        if(Integrands.size()<1) {
            IsZero = true;
            return;
        }
        
        if(Verbose > 1) cout << Color_HighLight << "  EpsEpExpands @ " << now() << RESET << endl;
        
        if(Verbose > 1) cout << "  \\--Collecting: " << Integrands.size() << " :> " << flush;
        map<ex, ex, ex_is_less> int_map;
        for(auto &item : Integrands) {
            if(item.is_zero()) continue;
            exset cts;
            item.find(CT(w), cts);
            if(cts.size() != 1) {
                cerr << "cts: " << cts << endl;
                cerr << "item: " << item << endl;
                throw Error("EpsEpExpands: CT size is NOT 1: ");
            }
            ex ct = (*(cts.begin())).subs(CT(w)==w);
            auto it = item.subs(CT(w)==1);
            int_map[it] = int_map[it]+ct;
        }
        
        Integrands.clear();
        Integrands.shrink_to_fit();
        for(auto kv : int_map) {
            if(kv.second.is_zero()) continue;
            Integrands.push_back(CT(kv.second) * kv.first);
        }
        if(Verbose > 1) cout << Integrands.size() << endl;

        auto res =
        GiNaC_Parallel(Integrands.size(), [&](int idx)->ex {
            // return { {two elements}, {two elements}, ...},
            // 1st: x-independent coefficient, expanded in ep/eps
            // 2nd: x-integrand
            auto item = Integrands[idx];
            if(item.is_zero()) return lst{ lst{0, 0} };
            exset cts;
            item.find(CT(w), cts);
            if(cts.size() != 1) {
                cerr << "cts: " << cts << endl;
                cerr << "item: " << item << endl;
                throw Error("EpsEpExpands: CT size is NOT 1: ");
            }
            ex ct = (*(cts.begin())).subs(CT(w)==w);
            auto it = item.subs(CT(w)==1);
            
            if(ct.has(epz)) {
                if(is_a<mul>(ct)) {
                    ex ct0 = 1, ct1 = 1;
                    for(auto ci : ct) {
                        if(ci.has(epz)) ct1 *= ci;
                        else ct0 *= ci;
                    }
                    ct = ct0;
                    it *= ct1;
                } else {
                    it *= ct;
                    ct = 1;
                }
            }
            
            if(it.has(epz)) it = mma_series(it,epz,0);
            auto cv_lst = mma_collect_lst(it, lst{epz, vs});
            if(cv_lst.nops()<1) return lst{ lst{ 0, 0} };
            
            lst para_res_lst;
            for(int i=0; i<cv_lst.nops();i++) {
                auto tmp = cv_lst.op(i).op(0);
                auto vc = cv_lst.op(i).op(1);
                //if(use_CCF) tmp = collect_common_factors(tmp);
                if(!tmp.has(eps) && !ct.has(eps)) {
                    if(tmp.has(epsID(w)) || ct.has(epsID(w))) {
                        throw Error("EpsEpExpands: epsID should be always multipled by eps!");
                    }
                    auto ct2 = vc * ct;
                    int ctN = epRank(ct2);
                    tmp = mma_series(tmp, ep, epN-ctN);
                    for(int di=tmp.ldegree(ep); (di<=tmp.degree(ep) && di<=epN-ctN); di++) {
                        auto intg = tmp.coeff(ep, di);
                        if(intg.has(ep)) {
                            throw Error("EpsEpExpands: ep found @ intg = " + ex2str(intg));
                        }
                        auto pref = mma_series(ct2, ep, epN-di);
                        if(pref.has(vs)) pref = mma_series(pref, vs, sN);
                        //if(use_CCF) intg = collect_common_factors(intg);
                        try{
                            intg = collect_common_factors(intg);
                        } catch(...) {
                            intg = collect_common_factors(intg.expand());
                        }
                        para_res_lst.append(lst{pref * pow(ep, di), intg});
                    }
                } else {
                    auto sct = vc * ct;
                    int sctN = epsRank(sct);
                    ex stmp = mma_series(tmp, eps, epsN-sctN);
                    for(int sdi=stmp.ldegree(eps); (sdi<=stmp.degree(eps) && sdi<=epsN-sctN); sdi++) {
                        tmp = stmp.coeff(eps, sdi);
                        //if(use_CCF) tmp = collect_common_factors(tmp);
                        if(tmp.has(eps)) {
                            throw Error("EpsEpExpands: eps found @ tmp = " + ex2str(tmp));
                        }
                        
                        auto cv_lst = mma_collect_lst(tmp, epsID(w));
                        auto ct2 = mma_series(sct, eps, epsN-sdi);
                        int ctN = epRank(ct2);
                        for(auto ti : cv_lst) { // Note: tmp is local
                            auto tmp = ti.op(0);
                            auto eps_ci = ti.op(1);
                            tmp = mma_series(tmp, ep, epN-ctN);
                            for(int di=tmp.ldegree(ep); (di<=tmp.degree(ep) && di<=epN-ctN); di++) {
                                auto intg = tmp.coeff(ep, di);
                                if(intg.has(ep)) {
                                    throw Error("EpsEpExpands: ep found @ intg = " + ex2str(intg));
                                }
                                auto pref = mma_series(ct2, ep, epN-di);
                                if(pref.has(vs)) pref = mma_series(pref, vs, sN);
                                //if(use_CCF) intg = collect_common_factors(intg);
                                try{
                                    intg = collect_common_factors(intg);
                                } catch(...) {
                                    intg = collect_common_factors(intg.expand());
                                }
                                para_res_lst.append(lst{eps_ci * pref * pow(eps, sdi) * pow(ep, di), intg});
                            }
                        }
                    }
                }
            }

            return para_res_lst;

        }, "EpsEp", !debug);
        
        if(Verbose > 1) cout << "  \\--Collecting: ";
        map<ex, ex, ex_is_less> int_pref;
        long long ncollect = 0;
        ex expr_nox = 0;
        for(auto &item : res) {
            ncollect += item.nops();
            for(auto &kv : ex_to<lst>(item)) {
                if(!kv.op(1).has(x(w))) expr_nox += kv.op(0) * kv.op(1);
                else int_pref[kv.op(1)] += kv.op(0);
            }
        }
        
        if(Verbose > 1) cout << ncollect << " :> " << flush;
        expResult.clear();
        expResult.shrink_to_fit();
        for(auto kv : int_pref) {
            if(kv.second.is_zero()) continue;
            expResult.push_back(lst{kv.second, kv.first});
        }
        if(!is_zero(expr_nox)) {
            expr_nox = expr_nox.subs(FTX(w1,w2)==1);
            expResult.push_back(lst{expr_nox, 1});
        }
        if(Verbose > 1) cout << expResult.size() << endl;
    }

}
