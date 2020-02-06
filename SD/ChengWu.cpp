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

vector<lst> SD::Binarize(lst fe, ex eqn) {
    auto xij = get_x_from(eqn);
    assert(xij.size()==2);
    ex xi = xij[0];
    ex xj = xij[1];
    ex ci = eqn.coeff(xi);
    ex cj = eqn.coeff(xj);
    assert((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0);
    
    vector<lst> ret;
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
    ret.push_back(fe1);
    
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
    ret.push_back(fe2);
    
    return ret;
}

bool SD::Partilize(ex f0, lst xs, lst &ret0, bool ext) {
    ex f = f0;
    for(auto xi : xs) {
        if(!f.has(xi) || f.degree(xi)!=1) continue;
        auto cxi = f.coeff(xi);
        if(!xPositive(cxi) && !xPositive(ex(0)-cxi)) continue;
        if(cxi.subs(x(wild())==1)<0) cxi = ex(0)-cxi;
        
        lst ret;
        ret.append(lst{xi, cxi});
        f = f.subs(xi==0);
        if(f.is_zero() || Partilize(f, xs, ret, ext)) {
            for(int i=0; i<ret.nops(); i++) ret0.append(ret.op(i));
            return true;
        }
        
        if(ext && (xPositive(f) || xPositive(ex(0)-f))) {
            ret0.append(lst{xi, cxi});
            if(f.subs(x(wild())==1)<0) f = ex(0)-f;
            ret0.append(lst{0, f});
            return true;
        }
    }
    return false;
}

void SD::ChengWu() {

    for(auto &fe : FunExp) {
        if(fe.nops()<3 || xPositive(fe.op(0).op(1)) || xPositive(0-fe.op(0).op(1))) continue;

        lst xs;
        for(int di=0; di<fe.op(2).nops(); di++) {
            auto delta = ex_to<lst>(fe.op(2).op(di));
            Projectivize(fe, delta);
            for(auto xi : delta) xs.append(xi);
        }
        
        auto ft = fe.op(0).op(1);
        ft = Factor(ft);
        while(true) {
            auto ft0 = ft;
            if(ft.match(pow(wild(1), wild(2)))) {
                ft = ft.op(0);
            } else if(is_a<mul>(ft)) {
                ex tmp = 1;
                for(auto fti : ft) {
                    if(xPositive(fti)) continue;
                    tmp = tmp * fti;
                }
                ft = tmp;
                if((ft-ft0).is_zero()) break;
                continue;
            }
            break;
        }

        lst ret;
        bool ok = Partilize(ft, xs, ret, false);
        if(!ok) {
            ret.remove_all();
            ok = Partilize(ft, xs, ret, true);
        }
        
        if(ok) {
            if(Verbose>10) cout << "  \\--" << WHITE << "Cheng-Wu @ size="  << ret.nops() << RESET << endl;
            auto ilast = ret.nops()-1;
            lst rm_xs;
            ex inv_det = 1;
            if(get_op(ret, ilast, 0).is_zero()) {
                ex xfi = x(x_free_index(fe));
                let_op(ret, ilast, 0, xfi);
                for(int i=ilast-1; i>=0; i--) {
                    let_op(ret, i, 1, get_op(ret,i,1)*xfi);
                }
                
                auto xs = get_op(fe,2,0,0);
                if(ChengWu_xsum) {
                    xs = 0;
                    for(auto xi : get_op(fe,2,0)) xs += xi;
                }
                let_op_append(fe, 2, 0, xfi);
                let_op_append(fe, 0, xs);
                let_op_append(fe, 0, xfi+xs);
                let_op_append(fe, 1, 1);
                let_op_append(fe, 1, -2);
                if(Verbose>10) cout << "    \\--" << WHITE << "Added xi1 = " << xfi << RESET << endl;
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
            lst x2y, num_xi_lst;
            for(auto ss : ret) {
                x2y.append(ss.op(0)==ss.op(1));
                if(!ss.op(1).has(x(wild()))) num_xi_lst.append(ss.op(0));
            }
                        
            lst re_xi_lst;
            for(int di=0; di<fe.op(2).nops(); di++) {
                auto delta = ex_to<lst>(fe.op(2).op(di));
                ex re_xi = 0;
                for(auto xi : delta) {
                    if(!rm_xs.has(xi) || num_xi_lst.has(xi)) {
                        re_xi = xi;
                        break;
                    }
                }
                if(re_xi.is_zero()) {
                    re_xi = x(x_free_index(fe));
                    auto xs = get_op(fe,2,di,0);
                    if(ChengWu_xsum) {
                        xs = 0;
                        for(auto xi : get_op(fe,2,di)) xs += xi;
                    }
                    let_op_append(fe, 2, di, re_xi);
                    let_op_append(fe, 0, xs);
                    let_op_append(fe, 0, re_xi+xs);
                    let_op_append(fe, 1, 1);
                    let_op_append(fe, 1, -2);
                    if(Verbose>10) cout << "    \\--" << WHITE << "Added xi2 = " << re_xi << RESET << endl;
                }
                re_xi_lst.append(re_xi);
            }

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

            for(int di=0; di<fe.op(2).nops(); di++) Projectivize(fe, fe.op(2).op(di), re_xi_lst.op(di));
        } else {
            // TODO: add more cases
        }
    }
    
    KillPowers();
}

}
