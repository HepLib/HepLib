#include "SD.h"
#include <cmath>

namespace HepLib {

void SD::Projectivize(ex &fe, ex delta, ex xsum) {
    symbol s;
    lst sRepl;
    for(int j=0; j<delta.nops(); j++) sRepl.append(delta.op(j)==delta.op(j)*s);
    if(xsum.is_zero()) {
        for(auto xi : delta) xsum += xi;
    }
    
    ex over_all_sn = 0;
    int nnn = fe.op(0).nops();
    for(int i=0; i<nnn; i++) {
        if(!fe.op(0).op(i).has(x(w))) continue;
        auto tmp = expand(fe.op(0).op(i));
        auto sn = tmp.subs(sRepl).degree(s);
        over_all_sn += sn*fe.op(1).op(i);
        lst items;
        if(is_a<add>(tmp)) {
            for(auto ii : tmp) items.append(ii);
        } else {
            items.append(tmp);
        }
        tmp = 0;
        for(auto ii : items) {
            auto sni = ii.subs(sRepl).degree(s);
            if(sni!=sn) tmp += ii * pow(xsum, sn-sni);
            else tmp += ii;
        }
        fe.let_op(0).let_op(i) = tmp;
    }

    over_all_sn = normal(over_all_sn+delta.nops());
    if(!over_all_sn.is_zero()) {
        let_op_append(fe, 0, xsum);
        let_op_append(fe, 1, ex(0)-over_all_sn);
    }
}

void SD::Scalelize(ex &fe, const ex xi, const ex cy) {
    if(is_a<lst>(xi)) Scalelize(fe, ex_to<lst>(xi), cy);
    else Scalelize(fe, lst{xi}, cy);
}
void SD::Scalelize(ex &fe, const lst xs, const ex cy) {
    lst x2y, y2x;
    for(auto xi : xs) {
        if(cy.has(xi)) {
            cerr << Color_Error << "Scalelize: cy has xi: " << cy << "/" << xi << RESET << endl;
            exit(1);
        }
        auto yi = xi.subs(x(w)==y(w));
        x2y.append(xi==cy*yi);
        y2x.append(yi==xi);
    }
    auto nnn = fe.op(0).nops();
    for(int i=0; i<nnn; i++) {
        if(!fe.op(0).op(i).has(x(w))) continue;
        auto tmp = fe.op(0).op(i).subs(x2y).subs(y2x);
        tmp = tmp.normal();
        tmp = tmp.numer_denom();
        if(tmp.op(1).subs(x(w)==1)<0) {
            tmp.let_op(0) = ex(0)-tmp.op(0);
            tmp.let_op(1) = ex(0)-tmp.op(1);
        }
        fe.let_op(0).let_op(i) = tmp.op(0);
        if(tmp.op(1)!=1) {
            let_op_append(fe, 0, tmp.op(1));
            let_op_append(fe, 1, ex(0)-fe.op(1).op(i));
        }
    }

    ex det = cy.normal();
    det = det.numer_denom();
    if(det.op(1).subs(x(w)==1)<0) {
        det.let_op(0) = ex(0)-det.op(0);
        det.let_op(1) = ex(0)-det.op(1);
    }
    auto xn = xs.nops();
    if(det.op(0)!=1) {
        let_op_append(fe, 0, det.op(0));
        let_op_append(fe, 1, xn);
    }
    if(det.op(1)!=1) {
        let_op_append(fe, 0, det.op(1));
        let_op_append(fe, 1, ex(0)-xn);
    }
}

vector<ex> SD::Binarize(ex const fe, ex const eqn) {
    vector<ex> add_to;
    auto xij = get_x_from(eqn);
    if(xij.size()!=2) {
        cerr << Color_Error << "Binarize: xij.size()!=2, " << xij << RESET << endl;
        exit(1);
    }
    ex xi = xij[0];
    ex xj = xij[1];
    ex ci = eqn.coeff(xi);
    ex cj = eqn.coeff(xj);
    if(!((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0)) {
        cerr << Color_Error << "Binarize: ((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0)" << RESET << endl;
        exit(1);
    }
    
    ci = abs(ci);
    cj = abs(cj);
    symbol yi,yj;
    // Part I: ci xi-cj xj>0, i.e., xi>cj/ci xj
    auto f1 = ex_to<lst>(fe.op(0));
    auto e1 = ex_to<lst>(fe.op(1));
    ex c1 = cj/ci;
    for(int i=0; i<f1.nops(); i++) {
        f1.let_op(i) = f1.op(i).subs(lst{xi==c1*yj/(1+c1)+yi,xj==yj/(1+c1)}).subs(lst{yi==xi,yj==xj});
    }
    if(e1.op(0)==1) f1.let_op(0) = f1.op(0)/(1+c1); // Jaccobi
    else if(e1.op(1)==1) f1.let_op(1) = f1.op(1)/(1+c1);
    else {
        f1.append(1/(1+c1));
        e1.append(1);
    }
    auto fe1 = fe;
    fe1.let_op(0) = f1;
    fe1.let_op(1) = e1;
    add_to.push_back(fe1);
    
    // Part II: ci xi-cj xj<0, i.e., i.e., xj>ci/cj xi
    auto f2 = ex_to<lst>(fe.op(0));
    auto e2 = ex_to<lst>(fe.op(1));
    ex c2 = ci/cj;
    for(int i=0; i<f2.nops(); i++) {
        f2.let_op(i) = f2.op(i).subs(lst{xj==c2*yi/(1+c2)+yj,xi==yi/(1+c2)}).subs(lst{yi==xi,yj==xj});
    }
    if(e2.op(0)==1) f2.let_op(0) = f2.op(0)/(1+c2); // Jaccobi
    else if(e2.op(1)==1) f2.let_op(1) = f2.op(1)/(1+c2);
    else {
        f2.append(1/(1+c2));
        e2.append(1);
    }
    auto fe2 = fe;
    fe2.let_op(0) = f2;
    fe2.let_op(1) = e2;
    add_to.push_back(fe2);
    return add_to;
}

/*
mode=0: x_i P_i, with P_i positive
mode=1: x_i P_i + G_0, with P_i and G_0 positive
mode=2: x_i P_i + G_0, with P-i positive, G_0 ~ (xm-xn)^n
mode=3: x_i P_i + x_0 G_0 + Q_0, with P-i positive, G_0 ~ (xm-xn)^n & Q_0 positive
mode=4: x_i P_i + x_0 G_0 + Q_0, with P-i positive, G_0/Q_0 ~ (xm-xn)^n
*/
bool SD::Partilize(ex f0, lst delta, lst &in_ret, int mode) {
    for(auto xi : delta) {
        ex f = f0;
        if(!f.has(xi) || f.degree(xi)!=1) continue;

        auto cxi = f.coeff(xi);
        f = f.subs(xi==0);
        int cxi_sgn = xSign(cxi);

        if(cxi_sgn!=0) {
            if(cxi_sgn<0) cxi = ex(0)-cxi;
        
            lst ret;
            ret.append(lst{xi, cxi});
            if(f.is_zero() || Partilize(f, delta, ret, mode)) { // mode=0
                for(int i=0; i<ret.nops(); i++) in_ret.append(ret.op(i));
                return true;
            }
            
            if((mode>0) && (xSign(f)!=0)) { // mode=1
                in_ret.append(lst{xi, cxi});
                if(f.subs(x(w)==1)<0) f = ex(0)-f;
                in_ret.append(lst{0, f});
                return true;
            } else if(mode>1) { // mode=2
                auto ff = Factor(f);
                lst fflst;
                if(is_a<mul>(ff)) {
                    for(auto item : ff) {
                        if(xSign(item)!=0) continue;
                        if(item.match(pow(w1, w2))) fflst.append(item.op(0));
                        else fflst.append(item);
                    }
                } else fflst.append(ff);
                
                if(fflst.nops()==1) {
                    symbol s;
                    ff = fflst.op(0).subs(x(w)==s*x(w));
                    if(get_x_from(ff).size()==2 && ff.degree(s)==1 && ff.ldegree(s)==1) {
                        in_ret.append(lst{xi, cxi});
                        in_ret.append(lst{0, fflst.op(0)});
                        return true;
                    }
                }
            }
        } else if(mode>2) {
            lst bilst;
            auto cc = Factor(cxi);
            lst cclst;
            if(is_a<mul>(cc)) {
                for(auto item : cc) {
                    if(xSign(item)!=0) continue;
                    if(item.match(pow(w1, w2))) cclst.append(item.op(0));
                    else cclst.append(item);
                }
            } else cclst.append(cc);
            if(cclst.nops()==1) {
                symbol s;
                cc = cclst.op(0).subs(x(w)==s*x(w));
                if(get_x_from(cc).size()==2 && cc.degree(s)==1 && cc.ldegree(s)==1) {
                    bilst.append(lst{0, cclst.op(0)});
                }
            } else continue;
            if(bilst.nops()!=1) continue;
            
            if(mode==3 && xSign(f)!=0) { // mode=3
                in_ret.append(bilst.op(0));
                return true;
            }
        
            if(mode==4) { // mode=4
                auto ff = Factor(f);
                lst fflst;
                if(is_a<mul>(ff)) {
                    for(auto item : ff) {
                        if(xSign(item)!=0) continue;
                        if(item.match(pow(w1, w2))) fflst.append(item.op(0));
                        else fflst.append(item);
                    }
                } else fflst.append(ff);
                if(fflst.nops()==1) {
                    symbol s;
                    ff = fflst.op(0).subs(x(w)==s*x(w));
                    if(get_x_from(ff).size()==2 && ff.degree(s)==1 && ff.ldegree(s)==1) {
                        bilst.append(lst{0, fflst.op(0)});
                    }
                } else continue;
                if(bilst.nops()!=2) continue;

                in_ret.append(bilst.op(0));
                bool ok = true;
                for(auto ii : get_x_from(bilst.op(0))) {
                    if(bilst.op(1).has(ii)) {
                        ok = false;
                        break;
                    }
                }
                if(ok) {
                    in_ret.append(bilst.op(1));
                    return true;
                }
            }
        }
        
        // TODO: other modes
    }
    return false;
}

void SD::ChengWu() {
    ChengWu(FunExp, Verbose);
}

// FunExp & Verbose are local
void SD::ChengWu(vector<ex> &FunExp, int Verbose) {
    vector<ex> FunExp2;
    for(auto fe : FunExp) {
        if(fe.nops()<3 || xSign(fe.op(0).op(1))!=0) {
            let_op_prepend(fe, 0, 1);
            let_op_prepend(fe, 1, 0);
            FunExp2.push_back(fe);
            continue;
        }
        auto deltas = fe.op(2);
        for(int di=0; di<deltas.nops(); di++) Projectivize(fe, deltas.op(di)); //make sure projective
        let_op_prepend(fe, 0, fe.op(0).op(1));
        let_op_prepend(fe, 1, 0);
        auto ret = ChengWu_Internal(fe, Verbose);
        for(auto item : ret) FunExp2.push_back(item);
    }
    
    // fe.op(0).op(0) : 1-ok, 2-nok
    // handle x_i P + Q, with Q: positive, P will apply Cheng-Wu 1st.
    FunExp.clear();
    for(auto fe : FunExp2) {
        if(is_zero(get_op(fe,0,0)-1)) {
            FunExp.push_back(fe);
            continue;
        }

        auto ft = RefinedFT(get_op(fe,0,2)).expand();
        auto xs = get_x_from(ft);

        for(auto xi : xs) {
            if(ft.degree(xi)==1 && xSign(ft.subs(xi==0))!=0) {
                auto fe2 = fe;
                let_op(fe2, 0, 0, ft.coeff(xi));
                auto ret1 = ChengWu_Internal(fe2, Verbose);
                for(auto item : ret1) {
                    if(get_op(item,0,0)!=1) goto inner_loop_end;
                }
                
                for(auto item : ret1) {
                    auto ft0 = get_op(item,0,2); // actual F-term
                    if(xSign(ft0)!=0) {
                        let_op(item, 0, 0, 1);
                        FunExp.push_back(item);
                        continue;
                    }
                    let_op(item, 0, 0, ft0); // actual F-term
                    if(Verbose>10) cout << "  \\--Cheng-Wu Subsection" << endl;
                    auto ret2 = ChengWu_Internal(item, Verbose);
                    for(auto item2 : ret2) FunExp.push_back(item2);
                }
                goto loop_end;
            }
            inner_loop_end: ;
        }
        FunExp.push_back(fe);
        loop_end: ;
    }
    
    // TODO: add more cases
    
    
    //remove the first item in op.(0) and op(1)
    for(auto &fe : FunExp) {
        let_op_remove_first(fe, 0);
        let_op_remove_first(fe, 1);
    }
}

// make sure ft in the first term ONLY appear in ChengWu.cpp
// input: ft = in_fe.op(0).op(0)
// ouput: in_op(0).op(0) replaced by 1-ok, 2-nok
vector<ex> SD::ChengWu_Internal(ex in_fe, int Verbose) {
    vector<ex> fe_lst, ret_lst;
    fe_lst.push_back(in_fe);
    while(true) {
        vector<ex> fe_lst2;
        for(int i=0; i<fe_lst.size(); i++) {
            auto fe = fe_lst[i];
            lst ret, delta;
            bool ok = false;
            
            auto ft = get_op(fe, 0, 0);
            // make sure, otherwise Projectivise may change things
            if(!get_op(fe, 1, 0).is_zero()) {
                cerr << Color_Error << "ChengWu_Internal: (!get_op(fe, 1, 0).is_zero())" << RESET << endl;
                exit(1);
            }
            ft = RefinedFT(ft);
            
            if(fe.nops()<3 || xSign(ft)!=0) {
                let_op(fe, 0, 0, 1);
                ret_lst.push_back(fe);
                ok = false; // so will not process following ok block
                goto ok_label;
            }
            
            for(int di=0; di<fe.op(2).nops(); di++) {
                delta = ex_to<lst>(fe.op(2).op(di));
                ok = Partilize(ft, delta, ret, 0);
                if(ok) {
                    if(Verbose>10) cout << "  \\--" << Color_HighLight << "Cheng-Wu @mode=0 and @size="  << ret.nops() << RESET << endl;
                    goto ok_label;
                }
                
                ret.remove_all();
                ok = Partilize(ft, delta, ret, 1);
                if(ok) {
                    if(Verbose>10) cout << "  \\--" << Color_HighLight << "Cheng-Wu @mode=1 and @size="  << ret.nops() << RESET << endl;
                    goto ok_label;
                }
                
                ret.remove_all();
                ok = Partilize(ft, delta, ret, 2);
                if(ok) {
                    auto bilst = Binarize(fe, get_op(ret, ret.nops()-1, 1));
                    for(auto item : bilst) fe_lst2.push_back(item);
                    ok = false;
                    if(Verbose>10) cout << "  \\--" << Color_HighLight << "Cheng-Wu @mode=2 and @size="  << ret.nops() << RESET << endl;
                    goto ok_label;
                }
                
                ret.remove_all();
                ok = Partilize(ft, delta, ret, 3);
                if(ok) {
                    auto bilst = Binarize(fe, get_op(ret, ret.nops()-1, 1));
                    for(auto item : bilst) fe_lst2.push_back(item);
                    ok = false;
                    if(Verbose>10) cout << "  \\--" << Color_HighLight << "Cheng-Wu @mode=3 and @size="  << ret.nops() << RESET << endl;
                    goto ok_label;
                }

                ret.remove_all();
                ok = Partilize(ft, delta, ret, 4);
                if(ok) {
                    auto eq1 = get_op(ret, ret.nops()-1, 1);
                    auto eq2 = get_op(ret, ret.nops()-2, 1);
                    for(auto item : Binarize(fe, eq1)) {
                        for(auto ii : Binarize(item, eq2)) fe_lst2.push_back(ii);
                    }
                    ok = false;
                    if(Verbose>10) cout << "  \\--" << Color_HighLight << "Cheng-Wu @mode=4 and @size="  << ret.nops() << RESET << endl;
                    goto ok_label;
                }
                
                //TODO: other cases
                
            }
            let_op(fe, 0, 0, 2);
            ret_lst.push_back(fe);
            
            ok_label:
            if(ok) {
                auto ilast = ret.nops()-1;
                lst rm_xs;
                ex inv_det = 1;
                if(get_op(ret,ilast,0).is_zero()) {
                    lst rep_xs;
                    for(int i=ilast-1; i>=0; i--) rep_xs.append(get_op(ret,i,0));
                    ex xfi=0;
                    for(auto xi : delta) {
                        if(!rep_xs.has(xi) && !get_op(ret,ilast,1).has(xi)) {
                            xfi = xi;
                            break;
                        }
                    }
                    if(is_zero(xfi)) {
                        ex xs0=0;
                        xfi = x(x_free_index(fe));
                        for(auto xi : delta) {
                            if(!rep_xs.has(xi)) {
                                xs0 = xi;
                                break;
                            }
                        }
                        if(xs0.is_zero()) {
                            cerr << Color_Error << "ChengWu_Internal: (xs0.is_zero())" << RESET << endl;
                            exit(1);
                        }
                        delta.append(xfi);
                        let_op_append(fe, 2, 0, xfi);
                        let_op_append(fe, 0, xs0);
                        let_op_append(fe, 0, xfi+xs0);
                        let_op_append(fe, 1, 1);
                        let_op_append(fe, 1, -2);
                    }
                    
                    let_op(ret, ilast, 0, xfi);
                    for(int i=ilast-1; i>=0; i--) let_op(ret, i, 1, get_op(ret,i,1)*xfi);
                }
                for(int i=ilast; i>=0; i--) {
                    auto xi = ret.op(i).op(0);
                    rm_xs.append(xi);
                    auto yi = xi.subs(x(w)==y(w));
                    auto s = ret.op(i).op(1);
                    inv_det *= s;
                    ret.let_op(i).let_op(1) = yi/s;
                    for(int j=i-1; j>=0; j--) {
                        ret.let_op(j) = ret.op(j).subs(xi==yi/s);
                    }
                }
                lst x2y;
                for(auto ss : ret) x2y.append(ss.op(0)==ss.op(1));
                
                auto nnn = fe.op(0).nops();
                for(int i=0; i<nnn; i++) {
                    if(!fe.op(0).op(i).has(x(w))) continue;
                    auto tmp = normal(fe.op(0).op(i).subs(x2y));
                    tmp = tmp.subs(y(w)==x(w));
                    auto num_den = numer_denom(tmp);
                    if(num_den.op(1).subs(x(w)==1)<0) {
                        num_den.let_op(0) = ex(0)-num_den.op(0);
                        num_den.let_op(1) = ex(0)-num_den.op(1);
                    }
                    fe.let_op(0).let_op(i) = num_den.op(0);
                    if(num_den.op(1)!=1) {
                        let_op_append(fe, 0, num_den.op(1));
                        let_op_append(fe, 1, ex(0)-fe.op(1).op(i));
                    }
                }

                inv_det = normal(inv_det.subs(y(w)==x(w)));
                auto idet_num_den = numer_denom(inv_det);
                if(idet_num_den.op(1).subs(x(w)==1)<0) {
                    idet_num_den.let_op(0) = ex(0)-idet_num_den.op(0);
                    idet_num_den.let_op(1) = ex(0)-idet_num_den.op(1);
                }
                if(idet_num_den.op(0)!=1) {
                    let_op_append(fe, 0, idet_num_den.op(0));
                    let_op_append(fe, 1, -1);
                }
                if(idet_num_den.op(1)!=1) {
                    let_op_append(fe, 0, idet_num_den.op(1));
                    let_op_append(fe, 1, 1);
                }
                
                ex re_xi = 0;
                for(auto xi : delta) {
                    if(!rm_xs.has(xi)) {
                        re_xi = xi;
                        break;
                    }
                }
                if(re_xi.is_zero()) {
                    if(rm_xs.nops()!=delta.nops()) {
                        cerr << Color_Error << "rm_xs = " << rm_xs << endl << "delta = " << delta << RESET << endl;
                        exit(1);
                    }
                    re_xi = rm_xs.op(0);
                }
                Projectivize(fe, delta, re_xi);
                
                auto new_ft = RefinedFT(get_op(fe, 0, 0));
                auto rxs = get_x_from(new_ft);
                lst xPos, xNeg;
                for(auto xi : rxs) {
                    if(new_ft.coeff(xi)>0) xPos.append(xi);
                    else xNeg.append(xi);
                }
                if(!is_zero(new_ft.subs(x(w)==1)-xPos.nops()+xNeg.nops()) && abs(new_ft)!=1) {
                    cerr << Color_Error << "ChengWu_Internal: new_ft = " << new_ft << RESET << endl;
                    exit(1);
                }
                
                if((xPos.nops()==1 && xNeg.nops()>1) || (xNeg.nops()==1 && xPos.nops()>1)) {
                    ex x0, x1, xsum;
                    if(xPos.nops()==1) {
                        x0 = xPos.op(0);
                        x1 = xNeg.op(0); // TODO: 任意性
                        xsum = ex(0)-new_ft.subs(x0==0);
                    } else {
                        x0 = xNeg.op(0);
                        x1 = xPos.op(0); // TODO: 任意性
                        xsum = new_ft.subs(x0==0);
                    }
                    Scalelize(fe, x0, xsum/x1);
                    auto bilst = Binarize(fe, x0-x1);
                    for(auto item : bilst) {
                        let_op(item, 0, 0, 1);
                        ret_lst.push_back(item);
                    }
                } else if(xPos.nops()>1 && xNeg.nops()>1) {
                    ex xPsum=0, xNsum=0;
                    for(auto xi : xPos) xPsum += xi;
                    for(auto xi : xNeg) xNsum += xi;
                    auto x0 = xPos.op(0); // TODO: 任意性
                    auto x1 = xNeg.op(0); // TODO: 任意性
                    Scalelize(fe, xPos, xNsum/x1);
                    Scalelize(fe, x1, xPsum/x0);
                    auto bilst = Binarize(fe, x0-x1);
                    for(auto item : bilst) {
                        let_op(item, 0, 0, 1);
                        ret_lst.push_back(item);
                    }
                } else if(xPos.nops()==1 && xNeg.nops()==1) {
                    auto bilst = Binarize(fe, new_ft);
                    for(auto item : bilst) {
                        let_op(item, 0, 0, 1);
                        ret_lst.push_back(item);
                    }
                } else {
                    let_op(fe, 0, 0, 1);
                    ret_lst.push_back(fe);
                }
            } // end of if(ok)
        }
        if(fe_lst2.size()<1) break;
        fe_lst = fe_lst2;
    }

    return ret_lst;
}

}
