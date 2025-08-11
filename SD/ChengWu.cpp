/**
 * @file
 * @brief Functions to use Cheng-Wu theorem to make F-term positive
 */

#include "SD.h"
#include <cmath>

namespace HepLib::SD {
    
    /**
     * @brief ChengWu for SecDec class
     */
    void SecDec::ChengWu(const ex & ft) {
        
        if(is_a<lst>(ft) && ft.nops()>1) {
            for(auto item : ft) {
                for(auto &fe : FunExp) {
                    let_op_prepend(fe, 0, item);
                    let_op_prepend(fe, 1, 0);
                }
            }
            
            int nt = ft.nops();
            for(int i=0; i<nt; i++) {
                auto FunExp2 = FunExp;
                FunExp.clear();
                FunExp.shrink_to_fit();
                for(auto fe : FunExp2) {
                    auto ft2 = fe.op(0).op(0);
                    if(!fe.op(1).op(0).is_equal(0)) throw Error("SecDec::ChengWu exponent is NOT 0.");
                    let_op_remove_first(fe, 0);
                    let_op_remove_first(fe, 1);
                    auto fes = ChengWu::Apply(fe, ft2);
                    for(auto fe2 : fes) {
                        let_op_remove_first(fe2, 0);
                        let_op_remove_first(fe2, 1);
                        FunExp.push_back(fe2);
                    }
                }
            }
        } else {
            if(is_a<lst>(ft)) FunExp = ChengWu::Apply(FunExp, ft.op(0));
            else FunExp = ChengWu::Apply(FunExp, ft);
            //remove the first item in op(0) and op(1)
            for(auto &fe : FunExp) {
                let_op_remove_first(fe, 0);
                let_op_remove_first(fe, 1);
            }
        }
    }
    
    /**
     * @brief ChengWu-rized on the vector of fe, 
     * note that 1st element of the output, 
     * which needs to be droped, its information is used to label the ChengWu is successful or NOT.
     * also check the Evaluate function for more information.
     * @param fe_vec the input vector of fe
     * @return ChengWu-rized vector of fe
     */
    exvector ChengWu::Apply(const exvector & fe_vec, const ex & ft) {
        exvector ret_fe_vec;
        for(auto fe : fe_vec) {
            ex cft = ft;
            if(is_zero(cft)) {
                if(fe.op(0).nops()>1) cft = fe.op(0).op(1);
                else cft = 1;
            }
            if(fe.nops()<3 || xSign(cft)!=0 || cft.has(WRA(w))) {
                let_op_prepend(fe, 0, 1);
                let_op_prepend(fe, 1, 0);
                ret_fe_vec.push_back(fe);
                continue;
            }
            auto deltas = fe.op(2);
            for(int di=0; di<deltas.nops(); di++) {
                Projectivize(fe, deltas.op(di)); // make sure projective
            }
            if(is_zero(ft)) cft = fe.op(0).op(1); // due to projective
            let_op_prepend(fe, 0, cft);
            let_op_prepend(fe, 1, 0);
            auto ret = Evaluate(fe);
            for(auto item : ret) ret_fe_vec.push_back(item);
        }

        return exvector(std::move(ret_fe_vec));
    }
    
    /**
     * @brief to check the input fe is Projective or NOT w.r.t. delta
     * @param fe is the { function list, exponet list }
     * @param delta is a list of x's in Delta function
     * @return true if fe is Projective
     */
    bool ChengWu::isProjective(const ex fe, const ex delta) {
        static symbol s;
        lst sRepl;
        for(int j=0; j<delta.nops(); j++) sRepl.append(delta.op(j)==delta.op(j)*s);
        ex over_all_sn = 0;
        int nnn = fe.op(0).nops();
        for(int i=0; i<nnn; i++) {
            auto func = fe.op(0).op(i);
            if(!func.has(x(w))) continue;
            func = expand_ex(func.subs(sRepl),s);
            auto sn = func.degree(s);
            if(sn!=func.ldegree(s)) return false;
            over_all_sn += sn*fe.op(1).op(i);
        }
        over_all_sn = normal(over_all_sn+delta.nops());
        return over_all_sn.is_zero();
    }

    /**
     * @brief to Projectivize the input fe, the fe will be updated
     * @param fe is the { function list, exponet list }
     * @param delta is a list of x's in Delta function
     * @param xsum_in is used to balance the powers, xsum refers to Delta(1-xsum), if xsum=0, then xsum=SUM(delta list)
     */
    void ChengWu::Projectivize(ex &fe, const ex delta, const ex xsum_in) {
        static symbol s;
        lst sRepl;
        for(int j=0; j<delta.nops(); j++) sRepl.append(delta.op(j)==delta.op(j)*s);
        ex xsum = xsum_in;
        if(xsum.is_zero()) {
            for(auto xi : delta) xsum += xi;
        }
        
        ex over_all_sn = 0;
        int nnn = fe.op(0).nops();
        for(int i=0; i<nnn; i++) {
            auto func = fe.op(0).op(i);
            if(!func.has(x(w))) continue;
            func = expand_ex(func.subs(sRepl),s);
            auto sn = func.degree(s);
            over_all_sn += sn*fe.op(1).op(i);
            if(!is_a<add>(func)) func = lst{func};
            ex tmp = 0;
            for(auto ii : func) {
                auto sni = ii.degree(s);
                if(sni!=sn) tmp += ii.subs(s==1) * pow(xsum, sn-sni);
                else tmp += ii.subs(s==1);
            }
            fe.let_op(0).let_op(i) = tmp;
        }

        over_all_sn = normal(over_all_sn+delta.nops());
        if(!over_all_sn.is_zero()) {
            let_op_append(fe, 0, xsum);
            let_op_append(fe, 1, ex(0)-over_all_sn);
        }
    }

    /**
     * @brief to Scalelize the input fe, the fe will be updated, make transformation xi -> cy * yi & yi->xi,  the determinant will be added to fe back
     * @param fe is the { function list, exponet list }
     * @param xi the variable
     * @param cy the coefficient
     */
    void ChengWu::Scalelize(ex &fe, const ex xi, const ex cy) {
        if(is_a<lst>(xi)) Scalelize(fe, ex_to<lst>(xi), cy);
        else Scalelize(fe, lst{xi}, cy);
    }
    
    /**
     * @brief to Scalelize the input fe, the fe will be updated
     * @param fe is the { function list, exponet list }
     * @param xs the list of { xi }, make transformation for each xi -> cy * yi & yi->xi,  the determinant will be added to fe back
     * @param cy the coefficient
     */
    void ChengWu::Scalelize(ex &fe, const lst xs, const ex cy) {
        lst x2y, y2x;
        for(auto xi : xs) {
            if(cy.has(xi)) {
                throw Error("Scalelize: cy has xi: " + ex2str(cy) + "/" + ex2str(xi));
            }
            auto yi = xi.subs(x(w)==y(w));
            x2y.append(xi==cy*yi);
            y2x.append(yi==xi);
        }
        auto nnn = fe.op(0).nops();
        for(int i=0; i<nnn; i++) {
            if(!fe.op(0).op(i).has(x(w))) continue;
            auto tmp = fe.op(0).op(i).subs(x2y).subs(y2x);
            tmp = exfactor(tmp);
            tmp = numer_denom(tmp);
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

        ex det = exfactor(cy);
        det = numer_denom(det);
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

    /**
     * @brief Binarize the input fe, w.r.t. the linear eqn
     * @param fe the { function list, exponet list }
     * @param eqn linear equation, line a xi + b xj, e.g. xi-xj can be devide into xi>xj and xi < xj, and then change xi/xj back to 0 to infinity
     * @return 2 elements of { function list, exponet list }, just like fe
     */
    exvector ChengWu::Binarize(ex const fe, ex const eqn) {
        exvector ovec;
        Binarize(fe, eqn, ovec);
        return exvector(std::move(ovec));
    }
    
    /**
     * @brief Binarize the input fe, w.r.t. the linear eqn
     * @param fe the { function list, exponet list }
     * @param eqn linear equation, line a xi + b xj, e.g. xi-xj can be devide into xi>xj and xi < xj, and then change xi/xj back to 0 to infinity
     * @param ovec will get updated
     */
    void ChengWu::Binarize(ex const fe, ex const eqn, exvector & ovec) {
        auto xij = get_x_from(eqn);
        if(xij.size()!=2) {
            throw Error("Binarize: xij.size()!=2, " + ex2str(xij));
        }
        ex xi = xij[0];
        ex xj = xij[1];
        ex ci = eqn.coeff(xi);
        ex cj = eqn.coeff(xj);
        if(!((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0)) {
            throw Error("Binarize: ((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0)");
        }
        
        bool pos1 = (ci>0);
        
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
        
        // return vector 
        if(pos1) {
            ovec.push_back(fe1);
            ovec.push_back(fe2);
        } else {
            ovec.push_back(fe2);
            ovec.push_back(fe1);
        }
    }
    
    /**
     * @brief Linearize w.r.t. F-term
     * @param ft the F-term
     * @param delta the delta list
     * @param xcs a list of format: { {x1,c1}, {x2,c2} ... }, will be updated by its append method
     * @return true for success
     */
    bool ChengWu::isLinearizable(const ex ft, const ex delta, lst & xcs) {
        for(auto xi : delta) {
            ex f = ft;
            if(!f.has(xi) || f.degree(xi)!=1) continue;
            auto cxi = f.coeff(xi);
            f = f.subs(xi==0);
            int cxi_sgn = xSign(cxi);

            if(cxi_sgn!=0) {
                if(cxi_sgn<0) cxi = ex(0)-cxi;
                lst ret;
                if(is_zero(f) || isLinearizable(f, delta, ret)) {
                    for(auto xc : ret) xcs.append(xc);
                    xcs.append(lst{xi, cxi});
                    return true;
                }
            } 
        }
        return false;
    }
    
    /**
     * @brief Linearize w.r.t. xcs_in
     * @param xcs_in from isLinerizable function
     * @param fe the input/ouput fe
     * @param ft other expression needs to be transformed
     */
    void ChengWu::Linearize(const lst xcs_in, ex & fe, ex & ft) {
        lst xcs = xcs_in;
        ex inv_det = 1;
        auto nnn = xcs.nops();
        for(int i=0; i<nnn; i++) {
            auto xi = xcs.op(i).op(0);
            auto s = xcs.op(i).op(1);
            auto yi = xi.subs(x(w)==y(w));
            inv_det *= s;
            xcs.let_op(i).let_op(1) = yi/s;
            for(int j=i+1; j<nnn; j++) {
                xcs.let_op(j) = xcs.op(j).subs(xi==yi/s);
            }
        }
        exmap x2y;
        for(auto ss : xcs) x2y[ss.op(0)]=ss.op(1);
        
        // rescaling ft
        if(ft.has(x(w))) ft = ft.subs(x2y).subs(y(w)==x(w));
        
        // rescaling each item
        nnn = fe.op(0).nops();
        for(int i=0; i<nnn; i++) {
            if(!fe.op(0).op(i).has(x(w))) continue;
            auto tmp = fe.op(0).op(i).subs(x2y);
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

        // determinant
        inv_det = inv_det.subs(y(w)==x(w));
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
    }

    /**
     * @brief isPartilize w.r.t. F-term, generized to isLinearizable
     * @param ft the F-term
     * @param delta the delta list
     * @param xcs a list of format: { {x1,c1}, {x2,c2} ... }, will be updated by its append method
     * @param mode 
     * in the following, the repeated i implies summation on it.
     * - mode=0: x_i P_i, with P_i same-sign
     * - mode=1: x_i P_i + G_0, with P_i and G_0 same-sign
     * - mode=2: x_i P_i + G_0, with P-i same-sign, G_0 ~ (xm-xn)^n
     * - mode=3: x_i P_i + x_0 G_0 + Q_0, with P-i same-sign, G_0 ~ (xm-xn)^n & Q_0 same-sign
     * - mode=4: x_i P_i + x_0 G_0 + Q_0, with P-i same-sign, G_0/Q_0 ~ (xm-xn)^n
     * @return 2 elements of { function list, exponet list }, just like fe
     */
    bool ChengWu::isPartilizable(const ex ft, const ex delta, lst &xcs, int mode) {
        for(auto xi : delta) {
            if(!ft.has(xi) || ft.degree(xi)!=1) continue;

            auto cxi = ft.coeff(xi);
            int cxi_sgn = xSign(cxi);
            auto re_ft = ft.subs(xi==0);

            if(cxi_sgn!=0) {
                if(cxi_sgn<0) cxi = ex(0)-cxi;
            
                lst ret;
                ret.append(lst{xi, cxi});
                if(is_zero(re_ft) || isPartilizable(re_ft, delta, ret, mode)) { // mode=0
                    for(int i=0; i<ret.nops(); i++) xcs.append(ret.op(i));
                    return true;
                }
                
                if((mode>0) && (xSign(re_ft)!=0)) { // mode=1: x_i P_i + G_0, with P_i and G_0 same-sign
                    xcs.append(lst{xi, cxi});
                    if(re_ft.subs(x(w)==1)<0) re_ft = ex(0)-re_ft;
                    xcs.append(lst{0, re_ft});
                    return true;
                } else if(mode>1) { // mode=2: x_i P_i + G_0, with P-i same-sign, G_0 ~ (xm-xn)^n
                    auto fflst = SecDec::XRefined_lst(re_ft);
                    if(fflst.nops()==1) {
                        symbol s;
                        auto ff = fflst.op(0).subs(x(w)==s*x(w));
                        if(get_x_from(ff).size()==2 && ff.degree(s)==1 && ff.ldegree(s)==1) {
                            xcs.append(lst{xi, cxi});
                            xcs.append(lst{0, fflst.op(0)});
                            return true;
                        }
                    }
                }
            } else if(mode>2) {
                lst bilst;
                auto cclst = SecDec::XRefined_lst(cxi);
                if(cclst.nops()==1) {
                    symbol s;
                    auto cc = cclst.op(0).subs(x(w)==s*x(w));
                    if(get_x_from(cc).size()==2 && cc.degree(s)==1 && cc.ldegree(s)==1) {
                        bilst.append(lst{0, cclst.op(0)});
                    }
                } else continue;
                if(bilst.nops()!=1) continue;
                
                if(mode==3 && xSign(re_ft)!=0) { 
                    // mode=3: x_i P_i + x_0 G_0 + Q_0, with P-i same-sign, G_0 ~ (xm-xn)^n & Q_0 same-sign
                    xcs.append(bilst.op(0));
                    return true;
                }
            
                if(mode==4) { 
                    // mode=4: x_i P_i + x_0 G_0 + Q_0, with P-i same-sign, G_0~(xm-xn)^n & Q_0~(xm'-xn')^n'
                    auto fflst = SecDec::XRefined_lst(re_ft);
                    if(fflst.nops()==1) {
                        symbol s;
                        auto ff = fflst.op(0).subs(x(w)==s*x(w));
                        if(get_x_from(ff).size()==2 && ff.degree(s)==1 && ff.ldegree(s)==1) {
                            bilst.append(lst{0, fflst.op(0)});
                        }
                    } else continue;
                    if(bilst.nops()!=2) continue;

                    bool ok = true;
                    for(auto ii : get_x_from(bilst.op(0))) {
                        if(bilst.op(1).has(ii)) {
                            ok = false;
                            break;
                        }
                    }
                    if(ok) {
                        xcs.append(bilst.op(0));
                        xcs.append(bilst.op(1));
                        return true;
                    }
                }
            }
            
            // TODO: other modes
        }
        return false;
    }
    
    /**
     * @brief Partilize function, generized to Linearize
     * @param xcs just from isPartilizable
     * @param delta_in the delta list
     * @param fe_in the original fe
     * @param ret_lst will be updated by push_back
     */
    void ChengWu::Partilize(const lst xcs, const lst delta_in, const ex fe_in, exvector & ret_lst) {
        if(Verbose>10) {
            cout << PRE << "\\--" << Color_HighLight << "ChengWu @xcs="  << xcs.nops() << RESET << endl;
        }
        lst ret = xcs;
        auto fe = fe_in;
        auto delta = delta_in;
        auto ilast = ret.nops()-1;
        ex inv_det = 1;
        if(is_zero(get_op(ret,ilast,0))) {
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
                for(auto xi : delta) {
                    if(!rep_xs.has(xi)) {
                        xs0 = xi;
                        break;
                    }
                }
                if(is_zero(xs0)) throw Error("Partilize: (xs0.is_zero())");
                xfi = x(x_free_index(fe));
                delta.append(xfi);
                int dlti = -1;
                for(int i=0; i<fe.op(2).nops(); i++) {
                    if(fe.op(2).op(i).has(xs0)) {
                        dlti = i;
                        break;
                    }
                }
                if(dlti<0) throw Error("Partilize: dlti<0 found.");
                let_op_append(fe, 2, dlti, xfi);
                let_op_append(fe, 0, xs0);
                let_op_append(fe, 0, xfi+xs0);
                let_op_append(fe, 1, 1);
                let_op_append(fe, 1, -2);
            }
            
            let_op(ret, ilast, 0, xfi);
            for(int i=ilast-1; i>=0; i--) let_op(ret, i, 1, get_op(ret,i,1)*xfi);
        }
        lst rm_xs;
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
            auto tmp = exfactor(fe.op(0).op(i).subs(x2y));
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

        inv_det = exfactor(inv_det.subs(y(w)==x(w)));
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
        
        if(!isProjective(fe,delta)) {
            ex re_xi = 0;
            for(auto xi : delta) {
                if(!rm_xs.has(xi)) {
                    re_xi = xi;
                    break;
                }
            }
            if(is_zero(re_xi)) {
                cout << "xcs = " << xcs << endl;
                cout << "fe = " << fe_in << endl;
                throw Error("Partilize: rm_xs = " + ex2str(rm_xs) + " & delta = " + ex2str(delta));
            }
            Projectivize(fe, delta, re_xi);
        }
        
        auto new_ft = SecDec::XRefined(get_op(fe, 0, 0));
        symbol ss;
        auto sft = new_ft.subs(x(w)==ss*x(w));
        if(sft.degree(ss)!=1 || sft.ldegree(ss)!=1) throw Error("Partilize: ldegree/degree is NOT 1.");
        auto rxs = get_x_from(new_ft);
        lst xPos, xNeg;
        ex sPos=0, sNeg=0;
        for(auto xi : rxs) {
            auto ci = new_ft.coeff(xi);
            if(ci>0) {
                xPos.append(xi); sPos += xi*ci;
            } else {
                xNeg.append(xi); sNeg -= xi*ci;
            }
        }
        if(is_zero(sPos*sNeg)) throw Error("Partilize: sPos * sNeg = 0");
        ex px = 0;
        ex nx = 0;
        for(int i=rm_xs.nops(); i>0; i--) {
            auto xi = rm_xs.op(i-1);
            if(is_zero(nx) && xNeg.has(xi)) nx = xi;
            else if(is_zero(px) && xPos.has(xi)) px = xi;
            if(!is_zero(px) && !is_zero(nx)) break;
        }
        if(sNeg<sPos) {
            ex s = sNeg/nx;
            if(s!=1) Scalelize(fe,xPos,s);
            s = sPos/px;
            if(s!=1) Scalelize(fe,nx,s);
        } else {
            ex s = sPos/px;
            if(s!=1) Scalelize(fe,xNeg,s);
            s = sNeg/nx;
            if(s!=1) Scalelize(fe,px,s);
        }
        ex eqn = px-nx;
        auto bilst = Binarize(fe, eqn);
        for(auto item : bilst) {
            if(!isProjective(item,delta)) throw Error("Partilize: something is wrong here.");
            let_op(item, 0, 0, 1);
            ret_lst.push_back(item);
        }
    }

    /**
     * @brief ChengWu Internal, make sure ft in the first term, ONLY appear in ChengWu.cpp
     * @param in_fe input fe, ft = in_fe.op(0).op(0), and in_fe.op(1).op(0) should be 0
     * @return vector of fe, and fe_op(0).op(0) replaced by 1-ok, 2-nok
     */
    exvector ChengWu::Evaluate(const ex & in_fe) {
        exvector fe_lst, ret_lst;
        fe_lst.push_back(in_fe);
        static int total_modes = 5;
        while(true) {
            exvector fe_lst2;
            for(int i=0; i<fe_lst.size(); i++) {
                auto fe = fe_lst[i];
                auto ft = get_op(fe, 0, 0);
                
                // make sure, otherwise Projectivise may change things
                if(!get_op(fe, 1, 0).is_zero()) {
                    throw Error("Evaluate: (!get_op(fe, 1, 0).is_zero())");
                }
                ft = SecDec::XRefined(ft);
                
                if(fe.nops()<3 || xSign(ft)!=0) {
                    let_op(fe, 0, 0, 1);
                    ret_lst.push_back(fe);
                    goto ok_label;
                }
                
                // loop modes and deltas
                for(int mode=0; mode<total_modes; mode++) {
                for(int di=0; di<fe.op(2).nops(); di++) {
                    lst delta = ex_to<lst>(fe.op(2).op(di));
                    if(delta.nops()<2) continue;
                    
                    lst ret;
                    if(mode==0 && isPartilizable(ft, delta, ret, mode)) {
                        Partilize(ret, delta, fe, ret_lst);
                        goto ok_label;
                    }
                    
                    ret.remove_all();
                    if(mode==1 && isPartilizable(ft, delta, ret, mode)) {
                        Partilize(ret, delta, fe, ret_lst);
                        goto ok_label;
                    }
                    
                    ret.remove_all();
                    if(mode==2 && isPartilizable(ft, delta, ret, mode)) {
                        auto bilst = Binarize(fe, get_op(ret, ret.nops()-1, 1));
                        for(auto item : bilst) fe_lst2.push_back(item);
                        goto ok_label;
                    }
                    
                    ret.remove_all();
                    if(mode==3 && isPartilizable(ft, delta, ret, mode)) {
                        auto bilst = Binarize(fe, get_op(ret, ret.nops()-1, 1));
                        for(auto item : bilst) fe_lst2.push_back(item);
                        goto ok_label;
                    }

                    ret.remove_all();
                    if(mode==4 && isPartilizable(ft, delta, ret, mode)) {
                        auto eq1 = get_op(ret, ret.nops()-1, 1);
                        auto eq2 = get_op(ret, ret.nops()-2, 1);
                        for(auto item : Binarize(fe, eq1)) {
                            for(auto ii : Binarize(item, eq2)) fe_lst2.push_back(ii);
                        }
                        goto ok_label;
                    }
                    
                    //TODO: other mode
                    
                }}
                let_op(fe, 0, 0, 2);
                ret_lst.push_back(fe);
                
                ok_label: ;
            }
            if(fe_lst2.size()<1) break;
            fe_lst = fe_lst2;
        }

        return exvector(std::move(ret_lst));
    }
    
    namespace {
        int maxn(const ex pols, const ex xi) {
            int max_xn = -1;
            for(auto item : pols) {
                auto cxn = item.degree(xi);
                if(max_xn<cxn) max_xn = cxn;
            }
            return max_xn;
        }
    }
    
    /**
     * @brief WickRotation, just check WRA exist or NOT to see successful or NOT. Still Experimental
     * @param fe_vec input fe vector
     * @return vector of fe
     */
    exvector ChengWu::WickRotation(const exvector & fe_vec) {
        exvector ret_vec, run_vec, run2_vec;
        run_vec = fe_vec;
        ReRun:
        for(auto fe : run_vec) {
            ex ft = fe.op(0).op(1);        
            auto xs = get_x_from(ft);
            
            // Case: x(cx!=0) + (c0!=0)
            for(auto xi : xs) {
                auto ftx = collect_ex(ft,xi);
                auto c0 = ftx.subs(xi==0);
                
                auto cx = ftx-c0;
                if(xSign(cx)!=0) {
                    if(xSign(c0)!=0) {
                        auto max_xn = maxn(fe.op(0),xi)+1;
                        auto wra = WRA(-xSign(cx) * Pi/max_xn);
                        fe.let_op(0) = fe.op(0).subs(xi==exp(I*wra)*xi);
                        if(fe.op(1).op(0)!=1) throw Error("WickRotation: fe.op(0).op(0)!=1.");
                        fe.let_op(0).let_op(0) = fe.op(0).op(0) * exp(I*wra);
                        ret_vec.push_back(fe);
                        goto next_fe;
                    } else {
                        auto item = SecDec::XRefined(c0);
                        auto xs2 = get_x_from(item);
                        for(auto xi2 : xs2) {
                            if(item.degree(xi2)!=1) continue;
                            auto xc0 = item.subs(xi2==0);
                            auto xc1 = item.coeff(xi2);
                            if(xSign(xc0)*xSign(xc1)<0) {
                                ex xx = 0;
                                for(auto xi3 : xs2) {
                                    if(xi3==xi || xi3==xi2) continue;
                                    xx = xi3;
                                    break;
                                } 
                                if(is_zero(xx)) continue;
                                Scalelize(fe,xi2,-xc0/xc1/xx);
                                auto bfe = Binarize(fe, xi2-xx);
                                for(auto item2 : bfe) run2_vec.push_back(item2);
                                goto next_fe;
                            }
                        }
                    }
                }
            }
            
            // Case: a*x^2+b*x+c with a*c<0 or b*c>0
            for(auto xi : xs) {
                auto ftx = collect_ex(ft,xi);
                if(ftx.degree(xi)!=2) continue;
                
                auto c0 = ftx.coeff(xi,0);
                auto c1 = ftx.coeff(xi,1);
                auto c2 = ftx.coeff(xi,2);
                
                if(xSign(c2) * xSign(c0) <0) {
                    if(xSign(c1)!=0) {
                        auto max_xn = maxn(fe.op(0),xi)+1;
                        auto wra = WRA(-xSign(c2) * Pi/max_xn);
                        fe.let_op(0) = fe.op(0).subs(xi==exp(I*wra)*xi);
                        fe.let_op(0).let_op(1) = c2*pow(xi,2)*exp(I*wra)+c1*xi+c0*exp(-I*wra);
                        if(fe.op(1).op(0)!=1) throw Error("WickRotation: fe.op(0).op(0)!=1.");
                        fe.let_op(0).let_op(0) = fe.op(0).op(0) * exp(I*wra*(1+fe.op(1).op(1)));
                        ret_vec.push_back(fe);
                        goto next_fe;
                    } else {
                        auto item = SecDec::XRefined(c1);
                        auto xs2 = get_x_from(item);
                        for(auto xi2 : xs2) {
                            if(item.degree(xi2)!=1) continue;
                            auto xc0 = item.subs(xi2==0);
                            auto xc1 = item.coeff(xi2);
                            if(xSign(xc0)*xSign(xc1)<0) {
                                ex xx = 0;
                                for(auto xi3 : xs2) {
                                    if(xi3==xi || xi3==xi2) continue;
                                    xx = xi3;
                                    break;
                                } 
                                if(is_zero(xx)) continue;
                                Scalelize(fe,xi2,-xc0/xc1/xx);
                                auto bfe = Binarize(fe, xi2-xx);
                                for(auto item2 : bfe) run2_vec.push_back(item2);
                                goto next_fe;
                            }
                        }
                    }
                }
                
                if(xSign(c1) * xSign(c0) >0) {
                    if(xSign(c2)!=0) {
                        auto max_xn = maxn(fe.op(0),xi)+1;
                        auto wra = WRA(xSign(c0) * Pi/max_xn);
                        fe.let_op(0) = fe.op(0).subs(xi==exp(I*wra)*xi);
                        fe.let_op(0).let_op(1) = c2*pow(xi,2)+c1*xi*exp(-I*wra)+c0*exp(-2*I*wra);
                        if(fe.op(1).op(0)!=1) throw Error("WickRotation: fe.op(0).op(0)!=1.");
                        fe.let_op(0).let_op(0) = fe.op(0).op(0) * exp(I*wra*(1+2*fe.op(1).op(1)));
                        ret_vec.push_back(fe);
                        goto next_fe;
                    } else {
                        auto item = SecDec::XRefined(c2);
                        auto xs2 = get_x_from(item);
                        for(auto xi2 : xs2) {
                            if(item.degree(xi2)!=1) continue;
                            auto xc0 = item.subs(xi2==0);
                            auto xc1 = item.coeff(xi2);
                            if(xSign(xc0)*xSign(xc1)<0) {
                                ex xx = 0;
                                for(auto xi3 : xs2) {
                                    if(xi3==xi || xi3==xi2) continue;
                                    xx = xi3;
                                    break;
                                } 
                                if(is_zero(xx)) continue;
                                Scalelize(fe,xi2,-xc0/xc1/xx);
                                auto bfe = Binarize(fe, xi2-xx);
                                for(auto item2 : bfe) run2_vec.push_back(item2);
                                goto next_fe;
                            }
                        }
                    }
                }
                
                auto c0x = Factor(c0);
                if(is_a<mul>(c0x)) {
                    ex cc = 1;
                    for(auto item : c0x) {
                        if(!is_a<numeric>(item) && !item.match(x(w))) cc *= item;
                    }
                    c0x = cc;
                }
                if(Factor(c0x).match(pow(w,2))) {
                    auto item = SecDec::XRefined(c0x);
                    auto xs2 = get_x_from(item);
                    for(auto xi2 : xs2) {
                        if(item.degree(xi2)!=1) continue;
                        auto xc0 = item.subs(xi2==0);
                        auto xc1 = item.coeff(xi2);
                        if(xSign(xc0)*xSign(xc1)<0) {
                            ex xx = 0;
                            for(auto xi3 : xs2) {
                                if(xi3==xi || xi3==xi2) continue;
                                xx = xi3;
                                break;
                            } 
                            if(is_zero(xx)) continue;
                            Scalelize(fe,xi2,-xc0/xc1/xx);
                            auto bfe = Binarize(fe, xi2-xx);
                            for(auto item2 : bfe) run2_vec.push_back(item2);
                            goto next_fe;
                        }
                    }
                }
            }
            
            ret_vec.push_back(fe);
            next_fe: ;
        }
        if(run2_vec.size()>0) {
            run_vec = run2_vec;
            run2_vec.clear();
            goto ReRun;
        }
        return exvector(std::move(ret_vec));
    }
    
}
