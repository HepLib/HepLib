#include "SD.h"
#include <cmath>

namespace HepLib {

void SD::Projectivize(lst &fe, lst delta, ex xsum) {
    symbol s;
    lst sRepl;
    for(int j=0; j<delta.nops(); j++) sRepl.append(delta.op(j)==delta.op(j)*s);
    if(xsum.is_zero()) {
        for(auto xi : delta) xsum += xi;
    }
    
    ex over_all_sn = 0;
    int nnn = fe.op(0).nops();
    for(int i=0; i<nnn; i++) {
        if(!fe.op(0).op(i).has(x(wild()))) continue;
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
void SD::Projectivize(lst &fe, ex delta, ex xsum) {
    Projectivize(fe, ex_to<lst>(delta), xsum);
}

void SD::Scalelize(lst &fe, ex xi, ex cyi) {
    auto nnn = fe.op(0).nops();
    auto yi = xi.subs(x(wild())==y(wild()));
    for(int i=0; i<nnn; i++) {
        if(!fe.op(0).op(i).has(x(wild()))) continue;
        auto tmp = fe.op(0).op(i).subs(xi==cyi*yi).subs(yi==xi);
        tmp = tmp.normal();
        tmp = tmp.numer_denom();
        if(tmp.op(1).subs(x(wild())==1)<0) {
            tmp.let_op(0) = ex(0)-tmp.op(0);
            tmp.let_op(1) = ex(0)-tmp.op(1);
        }
        fe.let_op(0).let_op(i) = tmp.op(0);
        if(tmp.op(1)!=1) {
            let_op_append(fe, 0, tmp.op(1));
            let_op_append(fe, 1, ex(0)-fe.op(1).op(i));
        }
    }
    ex det = cyi.normal();
    det = det.numer_denom();
    if(det.op(1).subs(x(wild())==1)<0) {
        det.let_op(0) = ex(0)-det.op(0);
        det.let_op(1) = ex(0)-det.op(1);
    }
    if(det.op(0)!=1) {
        let_op_append(fe, 0, det.op(0));
        let_op_append(fe, 1, 1);
    }
    if(det.op(1)!=1) {
        let_op_append(fe, 0, det.op(1));
        let_op_append(fe, 1, -1);
    }
}

vector<lst> SD::Binarize(lst const fe, ex const eqn) {
    vector<lst> add_to;
    auto xij = get_x_from(eqn);
    assert(xij.size()==2);
    ex xi = xij[0];
    ex xj = xij[1];
    ex ci = eqn.coeff(xi);
    ex cj = eqn.coeff(xj);
    assert((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0);
    
    ci = abs(ci);
    cj = abs(cj);
    symbol yi,yj;
    // Part I: ci xi-cj xj>0, i.e., xi>cj/ci xj
    auto f1 = fe.op(0);
    auto e1 = fe.op(1);
    ex c1 = cj/ci;
    for(int i=0; i<f1.nops(); i++) {
        f1.let_op(i) = f1.op(i).subs(lst{xi==c1*yj/(1+c1)+yi,xj==yj/(1+c1)}).subs(lst{yi==xi,yj==xj});
    }
    f1.let_op(0) = f1.op(0)/(1+c1); // Jaccobi
    auto fe1 = fe;
    fe1.let_op(0) = f1;
    fe1.let_op(1) = e1;
    add_to.push_back(fe1);
    
    // Part II: ci xi-cj xj<0, i.e., i.e., xj>ci/cj xi
    auto f2 = fe.op(0);
    auto e2 = fe.op(1);
    ex c2 = ci/cj;
    for(int i=0; i<f2.nops(); i++) {
        f2.let_op(i) = f2.op(i).subs(lst{xj==c2*yi/(1+c2)+yj,xi==yi/(1+c2)}).subs(lst{yi==xi,yj==xj});
    }
    f2.let_op(0) = f2.op(0)/(1+c2); // Jaccobi
    auto fe2 = fe;
    fe2.let_op(0) = f2;
    fe2.let_op(1) = e2;
    add_to.push_back(fe2);
    return add_to;
}

// Binarized to 2 items, one replace the input, one append to add_to
void SD::Binarize(lst &fe, ex const eqn, vector<lst> &add_to) {
    auto ret = Binarize(fe, eqn);
    auto fe1 = ret[0];
    auto fe2 = ret[1];
    fe.let_op(0) = fe1.op(0);
    fe.let_op(1) = fe1.op(1);
    add_to.push_back(fe2);
}

/*
mode=0: x_i P_i, with P_i positive
mode=1: x_i P_i + G_i, with P_i and G_i positive
mode=2: x_i P_i + G_i, with P-i positive, G_i ~ (xm-xn)^n
*/
bool SD::Partilize(ex f0, lst delta, lst &in_ret, int mode) {
    ex f = f0;
    for(auto xi : delta) {
        if(!f.has(xi) || f.degree(xi)!=1) continue;
        auto cxi = f.coeff(xi);
        int cxi_sgn = xSign(cxi);
        
        if(cxi_sgn!=0) {
            if(cxi_sgn<0) cxi = ex(0)-cxi;
        
            lst ret;
            ret.append(lst{xi, cxi});
            f = f.subs(xi==0);
            if(f.is_zero() || Partilize(f, delta, ret, mode)) {
                for(int i=0; i<ret.nops(); i++) in_ret.append(ret.op(i));
                return true;
            }
            
            if((mode>0) && (xSign(f)!=0)) {
                in_ret.append(lst{xi, cxi});
                if(f.subs(x(wild())==1)<0) f = ex(0)-f;
                in_ret.append(lst{0, f});
                return true;
            } else if(mode>1) {
                auto ff = Factor(f);
                lst bilst;
                if(is_a<mul>(ff)) {
                    for(auto item : ff) {
                        if(xSign(item)!=0) continue;
                        if(item.match(pow(wild(1), wild(2)))) bilst.append(item.op(0));
                        else bilst.append(item);
                    }
                } else {
                    bilst.append(ff);
                }
                
                if(bilst.nops()==1) {
                    symbol s;
                    ff = bilst.op(0).subs(x(wild())==s*x(wild()));
                    if(get_x_from(ff).size()==2 && ff.degree(s)==1 && ff.ldegree(s)==1) {
                        in_ret.append(lst{xi, cxi});
                        in_ret.append(lst{0, bilst.op(0)});
                        return true;
                    }
                }
            }
        }
        
        // TODO: other modes
    }
    return false;
}

void SD::ChengWu() {
    
    vector<lst> add_lst;
ChengWu_loop:
    for(auto item : add_lst) FunExp.push_back(item);
    add_lst.clear();
    add_lst.shrink_to_fit();
    for(auto &fe : FunExp) {
        if(fe.nops()<3 || xSign(fe.op(0).op(1))!=0) continue;
        
        for(int di=0; di<fe.op(2).nops(); di++) Projectivize(fe, fe.op(2).op(di));
        auto ft = fe.op(0).op(1);
        ft = Factor(ft);
        while(true) {
            auto ft0 = ft;
            if(ft.match(pow(wild(1), wild(2)))) {
                ft = ft.op(0);
            } else if(is_a<mul>(ft)) {
                ex tmp = 1;
                for(auto fti : ft) {
                    auto s = xSign(fti);
                    if(s>0) continue;
                    else if(s<0) tmp = ex(0)-tmp;
                    else tmp = tmp * fti;
                }
                ft = tmp;
                if((ft-ft0).is_zero()) break;
                continue;
            }
            break;
        }
        
        lst ret, delta;
        bool ok = false;
        for(int di=0; di<fe.op(2).nops(); di++) {
            delta = ex_to<lst>(fe.op(2).op(di));
            ok = Partilize(ft, delta, ret, 0);
            if(!ok) {
                ret.remove_all();
                ok = Partilize(ft, delta, ret, 1);
            } else {
                if(Verbose>10) cout << "  \\--" << WHITE << "Cheng-Wu @mode=0 and @size="  << ret.nops() << RESET << endl;
                goto ok_label;
            }

            if(!ok) {
                ret.remove_all();
                ok = Partilize(ft, delta, ret, 2);
                if(ok) {
                    Binarize(fe, get_op(ret, ret.nops()-1, 1), add_lst);
                    ok = false;
                    if(Verbose>10) cout << "  \\--" << WHITE << "Cheng-Wu @mode=2 and @size="  << ret.nops() << RESET << endl;
                    goto ok_label;
                }
            } else {
                if(Verbose>10) cout << "  \\--" << WHITE << "Cheng-Wu @mode=1 and @size="  << ret.nops() << RESET << endl;
                goto ok_label;
            }
        }
        
        ok_label:
        if(ok) {
            auto ilast = ret.nops()-1;
            lst rm_xs;
            ex inv_det = 1;
            if(get_op(ret, ilast, 0).is_zero()) {
                ex xfi = x(x_free_index(fe));
                let_op(ret, ilast, 0, xfi);
                for(int i=ilast-1; i>=0; i--) {
                    let_op(ret, i, 1, get_op(ret,i,1)*xfi);
                }
                auto xs0 = delta.op(0);
                if(ChengWu_xsum) {
                    xs0 = 0;
                    for(auto xi : delta) xs0 += xi;
                }
                delta.append(xfi);
                let_op_append(fe, 2, 0, xfi);
                let_op_append(fe, 0, xs0);
                let_op_append(fe, 0, xfi+xs0);
                let_op_append(fe, 1, 1);
                let_op_append(fe, 1, -2);
                if(Verbose>10) cout << "    \\--" << WHITE << "Added " << xfi << " to " << delta << RESET << endl;
            }
            for(int i=ilast; i>=0; i--) {
                auto xi = ret.op(i).op(0);
                rm_xs.append(xi);
                auto yi = xi.subs(x(wild())==y(wild()));
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
                if(!fe.op(0).op(i).has(x(wild()))) continue;
                auto tmp = normal(fe.op(0).op(i).subs(x2y));
                tmp = tmp.subs(y(wild())==x(wild()));
                auto num_den = numer_denom(tmp);
                if(num_den.op(1).subs(x(wild())==1)<0) {
                    num_den.let_op(0) = ex(0)-num_den.op(0);
                    num_den.let_op(1) = ex(0)-num_den.op(1);
                }
                fe.let_op(0).let_op(i) = num_den.op(0);
                if(num_den.op(1)!=1) {
                    let_op_append(fe, 0, num_den.op(1));
                    let_op_append(fe, 1, ex(0)-fe.op(1).op(i));
                }
            }

            inv_det = normal(inv_det.subs(y(wild())==x(wild())));
            auto idet_num_den = numer_denom(inv_det);
            if(idet_num_den.op(1).subs(x(wild())==1)<0) {
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
                re_xi = rm_xs.op(0);
                cout << "ret = " << ret << endl;
                cout << RED <<  "use re_xi = " << re_xi << RESET << endl;
            }
            Projectivize(fe, delta, re_xi);
        } else {
            // TODO: add more cases
        }
    }
    if(add_lst.size()>0) goto ChengWu_loop;
    
    KillPowers();
}

}
