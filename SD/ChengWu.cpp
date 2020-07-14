 /**
 * @file
 * @brief Functions to use Cheng-Wu theorem to make F-term positive
 * @author F. Feng
 * @version 1.0.0
 * @date 2020-04-21
 */

#include "SD.h"
#include <cmath>

namespace HepLib::SD {

    /**
     * @brief ChengWu for SecDec class
     * @param sub_cw     
     */
    void SecDec::ChengWu(bool sub_cw) {
        ChengWu::Apply(FunExp, sub_cw);
    }
    
    /**
     * @brief ChengWu, note that FunExp is now local
     * @param FunExp will be updated
     * @param sub_cw     
     */
    void ChengWu::Apply(vector<ex> &FunExp, bool sub_cw) {
        vector<ex> FunExp2;
        for(auto fe : FunExp) {
            if(fe.nops()<3 || xSign(fe.op(0).op(1))!=0 || fe.op(0).op(1).has(WRA(w))) {
                let_op_prepend(fe, 0, 1);
                let_op_prepend(fe, 1, 0);
                FunExp2.push_back(fe);
                continue;
            }
            auto deltas = fe.op(2);
            for(int di=0; di<deltas.nops(); di++) Projectivize(fe, deltas.op(di)); //make sure projective
            let_op_prepend(fe, 0, fe.op(0).op(1));
            let_op_prepend(fe, 1, 0);
            auto ret = Evaluate(fe);
            for(auto item : ret) FunExp2.push_back(item);
        }
        
        // fe.op(0).op(0) : 1-ok, 2-nok
        // handle x_i P + Q, with Q: positive, P will apply Cheng-Wu 1st.
        FunExp.clear();
        for(auto fe : FunExp2) {
            if(!sub_cw || is_zero(get_op(fe,0,0)-1)) {
                FunExp.push_back(fe);
                continue;
            }

            auto ft = SecDec::RefinedFT(get_op(fe,0,2)).expand();
            auto xs = get_x_from(ft);

            for(auto xi : xs) {
                if(ft.degree(xi)==1 && xSign(ft.subs(xi==0))!=0) {
                    auto fe2 = fe;
                    let_op(fe2, 0, 0, ft.coeff(xi));
                    auto ret1 = Evaluate(fe2);
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
                        auto ret2 = Evaluate(item);
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

    /**
     * @brief to Projectivize the input fe, the fe will be updated
     * @param fe is the { function list, exponet list }
     * @param delta is a list of x's in Delta function
     * @param xsum is used to balance the powers, xsum refers to Delta(1-xsum), if xsum=0, then xsum=SUM(delta list)
     */
    void ChengWu::Projectivize(ex &fe, ex delta, ex xsum) {
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
                cerr << ErrColor << "Scalelize: cy has xi: " << cy << "/" << xi << RESET << endl;
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

    /**
     * @brief Binarize the input fe, w.r.t. the linear eqn
     * @param fe the { function list, exponet list }
     * @param eqn linear equation, line a xi + b xj, e.g. xi-xj can be devide into xi>xj and xi < xj, and then change xi/xj back to 0 to infinity
     * @return 2 elements of { function list, exponet list }, just like fe
     */
    vector<ex> ChengWu::Binarize(ex const fe, ex const eqn) {
        vector<ex> add_to;
        auto xij = get_x_from(eqn);
        if(xij.size()!=2) {
            cerr << ErrColor << "Binarize: xij.size()!=2, " << xij << RESET << endl;
            exit(1);
        }
        ex xi = xij[0];
        ex xj = xij[1];
        ex ci = eqn.coeff(xi);
        ex cj = eqn.coeff(xj);
        if(!((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0)) {
            cerr << ErrColor << "Binarize: ((ci*xi+cj*xj-eqn).is_zero() && is_a<numeric>(ci * cj) && (ci*cj)<0)" << RESET << endl;
            exit(1);
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
            add_to.push_back(fe1);
            add_to.push_back(fe2);
        } else {
            add_to.push_back(fe2);
            add_to.push_back(fe1);
        }
        return add_to;
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
     * - mode=0: x_i P_i, with P_i positive
     * - mode=1: x_i P_i + G_0, with P_i and G_0 positive
     * - mode=2: x_i P_i + G_0, with P-i positive, G_0 ~ (xm-xn)^n
     * - mode=3: x_i P_i + x_0 G_0 + Q_0, with P-i positive, G_0 ~ (xm-xn)^n & Q_0 positive
     * - mode=4: x_i P_i + x_0 G_0 + Q_0, with P-i positive, G_0/Q_0 ~ (xm-xn)^n
     * @return 2 elements of { function list, exponet list }, just like fe
     */
    bool ChengWu::isPartilizable(const ex ft, const ex delta, lst &xcs, int mode) {
        for(auto xi : delta) {
            ex f = ft;
            if(!f.has(xi) || f.degree(xi)!=1) continue;

            auto cxi = f.coeff(xi);
            f = f.subs(xi==0);
            int cxi_sgn = xSign(cxi);

            if(cxi_sgn!=0) {
                if(cxi_sgn<0) cxi = ex(0)-cxi;
            
                lst ret;
                ret.append(lst{xi, cxi});
                if(f.is_zero() || isPartilizable(f, delta, ret, mode)) { // mode=0
                    for(int i=0; i<ret.nops(); i++) xcs.append(ret.op(i));
                    return true;
                }
                
                if((mode>0) && (xSign(f)!=0)) { // mode=1: x_i P_i + G_0, with P_i and G_0 positive
                    xcs.append(lst{xi, cxi});
                    if(f.subs(x(w)==1)<0) f = ex(0)-f;
                    xcs.append(lst{0, f});
                    return true;
                } else if(mode>1) { // mode=2: x_i P_i + G_0, with P-i positive, G_0 ~ (xm-xn)^n
                    auto fflst = SecDec::RefinedFT_lst(f);
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
                auto cclst = SecDec::RefinedFT_lst(cxi);
                if(cclst.nops()==1) {
                    symbol s;
                    auto cc = cclst.op(0).subs(x(w)==s*x(w));
                    if(get_x_from(cc).size()==2 && cc.degree(s)==1 && cc.ldegree(s)==1) {
                        bilst.append(lst{0, cclst.op(0)});
                    }
                } else continue;
                if(bilst.nops()!=1) continue;
                
                if(mode==3 && xSign(f)!=0) { 
                    // mode=3: x_i P_i + x_0 G_0 + Q_0, with P-i positive, G_0 ~ (xm-xn)^n & Q_0 positive
                    xcs.append(bilst.op(0));
                    return true;
                }
            
                if(mode==4) { 
                    // mode=4: x_i P_i + x_0 G_0 + Q_0, with P-i positive, G_0/Q_0 ~ (xm-xn)^n
                    auto fflst = SecDec::RefinedFT_lst(f);
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
            cout << "  \\--" << Color_HighLight << "Cheng-Wu @xcs="  << xcs.nops() << RESET << endl;
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
                xfi = x(x_free_index(fe));
                for(auto xi : delta) {
                    if(!rep_xs.has(xi)) {
                        xs0 = xi;
                        break;
                    }
                }
                if(xs0.is_zero()) throw Error("Partilize: (xs0.is_zero())");
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
        if(is_zero(re_xi)) throw Error("Partilize: rm_xs = " + ex2str(rm_xs) + "delta = " + ex2str(delta));
        Projectivize(fe, delta, re_xi);
        
        auto new_ft = SecDec::RefinedFT(get_op(fe, 0, 0));
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
        ex s = sNeg/xNeg.op(0);
        if(s!=1) Scalelize(fe,xPos,s);
        s = sPos/xPos.op(0);
        if(s!=1) Scalelize(fe,xNeg.op(0),s);
        ex eqn = xPos.op(0)-xNeg.op(0);
        auto bilst = Binarize(fe, eqn);
        for(auto item : bilst) {
            let_op(item, 0, 0, 1);
            ret_lst.push_back(item);
        }
    }

    /**
     * @brief ChengWu Internal, make sure ft in the first term, ONLY appear in ChengWu.cpp
     * @param in_fe input fe, ft = in_fe.op(0).op(0), and in_fe.op(1).op(0) should be 0
     * @return vector of fe, and fe_op(0).op(0) replaced by 1-ok, 2-nok
     */
    vector<ex> ChengWu::Evaluate(ex in_fe) {
        vector<ex> fe_lst, ret_lst;
        fe_lst.push_back(in_fe);
        static int total_modes = 5;
        while(true) {
            vector<ex> fe_lst2;
            for(int i=0; i<fe_lst.size(); i++) {
                auto fe = fe_lst[i];
                auto ft = get_op(fe, 0, 0);
                
                // make sure, otherwise Projectivise may change things
                if(!get_op(fe, 1, 0).is_zero()) {
                    throw Error("Evaluate: (!get_op(fe, 1, 0).is_zero())");
                }
                ft = SecDec::RefinedFT(ft);
                
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
                        Partilize(ret, delta, fe, fe_lst);
                        goto ok_label;
                    }
                    
                    ret.remove_all();
                    if(mode==1 && isPartilizable(ft, delta, ret, mode)) {
                        Partilize(ret, delta, fe, fe_lst);
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

        return ret_lst;
    }

}
