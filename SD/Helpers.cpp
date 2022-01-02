/**
 * @file
 * @brief Helpers Functions for SD
 */
 
#include "SD.h"

namespace HepLib::SD {

    /*-----------------------------------------------------*/
    // Global Functions
    /*-----------------------------------------------------*/

    exvector get_xy_from(ex pol) {
        exset xyset;
        bool ok = pol.find(x(w), xyset);
        if(!ok) {
            ok = pol.find(y(w), xyset);
            if(!ok) return exvector();
        }
        exvector xys(xyset.size());
        copy(xyset.begin(), xyset.end(), xys.begin());
        sort(xys.begin(), xys.end(), [](const auto &a, const auto &b){
            return normal((b-a)).subs(lst{ x(w)==w, y(w)==w }).info(info_flags::positive);
        });
        return exvector(std::move(xys));
    }

    exvector get_x_from(ex pol) {
        exset xset;
        bool ok = pol.find(x(w), xset);
        if(!ok) {
            exvector xs(0);
            return xs;
        }
        exvector xs(xset.size());
        copy(xset.begin(), xset.end(), xs.begin());
        sort(xs.begin(), xs.end(), [](const auto &a, const auto &b){
            return normal((b-a)).subs(lst{ x(w)==w }).info(info_flags::positive);
        });
        return exvector(std::move(xs));
    }

    exvector get_y_from(ex pol) {
        exset yset;
        bool ok = pol.find(y(w), yset);
        if(!ok) {
            exvector ys(0);
            return ys;
        }
        exvector ys(yset.size());
        copy(yset.begin(), yset.end(), ys.begin());
        sort(ys.begin(), ys.end(), [](const auto &a, const auto &b){
            return normal((b-a)).subs(lst{ y(w)==w }).info(info_flags::positive);
        });
        return exvector(std::move(ys));
    }
    
    exvector get_z_from(ex pol) {
        exset zset;
        bool ok = pol.find(z(w), zset);
        if(!ok) {
            exvector zs(0);
            return zs;
        }
        exvector zs(zset.size());
        copy(zset.begin(), zset.end(), zs.begin());
        sort(zs.begin(), zs.end(), [](const auto &a, const auto &b){
            return normal((b-a)).subs(lst{ z(w)==w }).info(info_flags::positive);
        });
        return exvector(std::move(zs));
    }

    exvector get_pl_from(ex pol) {
        exset plset;
        bool ok = pol.find(PL(w), plset);
        if(!ok) {
            exvector pls(0);
            return pls;
        }
        exvector pls(plset.size());
        copy(plset.begin(), plset.end(), pls.begin());
        sort(pls.begin(), pls.end(), [](const auto &a, const auto &b){
            return normal((b-a)).subs(lst{ PL(w)==w }).info(info_flags::positive);
        });
        return exvector(std::move(pls));
    }

    /*-----------------------------------------------------*/
    // VE
    /*-----------------------------------------------------*/
    ex VESimplify(ex expr, int epN, int epsN) {
        auto expr1 = EvalF(expr);
        if(expr1.has(eps) && !expr1.is_polynomial(eps)) expr1 = series_ex(expr1, eps, epsN);
        if(expr1.has(ep) && !expr1.is_polynomial(ep)) expr1 = series_ex(expr1, ep, epN);
        expr1 = collect_ex(expr1, lst{eps,ep});
        ex ret = 0;
        for(int si=expr1.ldegree(eps); si<=epsN; si++) {
        for(int i=expr1.ldegree(ep); i<=epN; i++) {
        
            auto ccRes = expr1.coeff(eps,si).coeff(ep,i);
            ccRes = NN(ccRes).expand();
            if(!is_a<add>(ccRes)) ccRes = lst{ ccRes };
            
            exmap pvmap;
            for(auto item : ccRes) {
                ex cc = 1, cf = 1;
                if(is_a<mul>(item)) {
                    cc = 1; 
                    cf = 1;
                    for(auto ii : item) {
                        if(is_a<numeric>(ii) || ii.match(VE(w1, w2))) {
                            cc *= ii;
                        } else {
                            cf *= ii;
                        }
                    }
                } else if(is_a<numeric>(item) || item.match(VE(w1, w2))) {
                    cc = item;
                } else {
                    cf = item;
                }
                pvmap[cf] += cc;
            }
        
            for(auto kv : pvmap) {
                auto vf = kv.first;
                auto cvs = collect_lst(kv.second, VE(w1,w2));
                ex vIR=0, eI2 = 0, eR2 = 0;
                for(auto cv : cvs) {
                    auto co = NN(cv.op(0));
                    auto ve = cv.op(1);
                    if(!is_a<numeric>(co)) {
                        cout << cv << endl;
                        throw Error("VESimplify: coefficient of VE is not numeric.");
                    }
                    if(abs(co)>numeric("1E300")) return NaN;
                    vIR += co * ve.op(0);
                    numeric nco = ex_to<numeric>(co);
                    ex ee = ve.op(1) * ve.op(1);
                    eR2 += nco.real_part() * nco.real_part() * ee;
                    eI2 += nco.imag_part() * nco.imag_part() * ee;
                }
                if(!is_zero(eR2) || !is_zero(vIR))
                    ret += VE(ex_to<numeric>(vIR).real_part(), sqrt(eR2)) * pow(eps,si) * pow(ep,i) * vf;
                if(!is_zero(eI2) || !is_zero(vIR))
                    ret += VE(ex_to<numeric>(vIR).imag_part(), sqrt(eI2)) * pow(eps,si) * pow(ep,i) * vf * I;
            }
        }}
        
        return ret.collect(lst{eps,ep}, true);
    }

    ex VEResult(ex expr) {
        return expr.subs(VE(0,0)==0).subs(VE(w1,w2)==VEO(w1,w2));
    }
    
    ex VEResult2(ex expr) {
        return expr.subs(VE(0,0)==0).subs(VE(w1,w2)==VEO2(w1,w2));
    }
    
    ex VEMaxErr(ex expr) {
        
        auto ccRes = expr.expand();
        lst ccResList;
        if(is_a<add>(ccRes)) {
            for(auto item : ccRes) ccResList.append(item);
        } else {
            ccResList.append(ccRes);
        }
            
        ccRes = -100;
        for(auto item : ccResList) {
            ex ntmp;
            if(is_a<mul>(item)) {
                ntmp = 1;
                for(auto ii : item) {
                    if(is_a<numeric>(ii) || ii.match(VE(w1, w2))) {
                        ntmp *= ii;
                    } 
                }
            } else if(is_a<numeric>(item) || item.match(VE(w1, w2))) {
                ntmp = item;
            }
            ntmp = NN(abs(ntmp.subs(VE(w1,w2)==w2)));
            if(ntmp>ccRes) ccRes = ntmp;
        }
        
        return ccRes;
    }


    /*-----------------------------------------------------*/
    // Heplers used in SD
    /*-----------------------------------------------------*/
    int x_free_index(ex expr) {
        auto xs = get_x_from(expr);
        for(int i=0; i<=xs.size(); i++) {
            bool ok = true;
            for(auto xi : xs) {
                if(xi.is_equal(x(i))) {
                    ok = false;
                    break;
                }
            }
            if(ok) return i;
        }
        return -1;
    }

    int y_free_index(ex expr) {
        auto ys = get_y_from(expr);
        for(int i=0; i<=ys.size(); i++) {
            bool ok = true;
            for(auto yi : ys) {
                if(yi.is_equal(y(i))) {
                    ok = false;
                    break;
                }
            }
            if(ok) return i;
        }
        return -1;
    }

    int epRank(ex expr_in) {
        if(!expr_in.has(ep)) return 0;
        int p = -5;
        auto expr = collect_ex(expr_in, ep);
        while(true) {
            auto tmp = series_to_poly(expr.series(ep, p));
            if(!tmp.is_zero()) {
                tmp = collect_ex(tmp, ep);
                return tmp.ldegree(ep);
            } else p++;
        }
        throw Error("epRank error!");
    }

    int epsRank(ex expr_in) {
        if(!expr_in.has(eps)) return 0;
        int p = -5;
        auto expr = collect_ex(expr_in, eps);
        while(true) {
            auto tmp = series_to_poly(expr.series(eps, p));
            if(!tmp.is_zero()) {
                tmp = collect_ex(tmp, eps);
                return tmp.ldegree(eps);
            } else p++;
        }
        throw Error("epsRank error!");
    }

    int vsRank(ex expr_in) {
        if(!expr_in.has(vs)) return 0;
        int p = -5;
        auto expr = collect_ex(expr_in, vs);
        while(true) {
            auto tmp = series_to_poly(expr.series(vs, p));
            if(!tmp.is_zero()) {
                tmp = collect_ex(tmp, vs);
                return tmp.ldegree(vs);
            } else p++;
        }
        throw Error("vsRank error!");
    }

/*-----------------------------------------------------*/
// IntegratorBase::inDQMP
/*-----------------------------------------------------*/
int IntegratorBase::inDQMP(qREAL const *x) {
    unsigned int xdim = XDim;
    
    if(xdim<=MPXDim) return 3;
    qREAL xmin = 100;
    for(int i=0; i<xdim; i++) {
        if(x[i] < xmin) xmin = x[i];
    }
    if(xmin < MPXLimit) return 3;
    
    if(FT!=NULL) {
        qREAL ft = 1E50;
        static FT_Type last_ft = NULL;
        static qREAL ft0;
        if(last_ft!=FT) {
            qREAL x0[xdim];
            #ifdef _USE_FLOAT128
            for(int i=0; i<xdim; i++) x0[i]=0.521Q;
            #else
            for(int i=0; i<xdim; i++) x0[i]=0.521L;
            #endif

            #ifdef _USE_FLOAT128
            qREAL ft0 = fabsq(FT(x0, Parameter));
            #else
            qREAL ft0 = fabsl(FT(x0, Parameter));
            #endif
            if(ft0<1E-50) ft0 = 1;
            last_ft = FT;
        } 
        #ifdef _USE_FLOAT128
        ft = fabsq(FT(x, Parameter));
        #else
        ft = fabsl(FT(x, Parameter));
        #endif
        ft = ft/ft0;
        if(ft<MPFLimit) return 3;
        else if(ft<QFLimit) return 2;
    }
    
    if(xdim <= QXDim || xmin < QXLimit) return 2;
    
    return 1;
}

ex SecDec::VEResult() {
    return HepLib::SD::VEResult(ResultError);
}

void SecDec::VEPrint(bool endlQ) {
    ex expr = HepLib::SD::VEResult(ResultError);
    for(int i=expr.ldegree(eps); i<=expr.degree(eps); i++) {
        ex exp1 = expr.coeff(eps, i);
        for(int j=expr.ldegree(ep); j<=expr.degree(ep); j++) {
            cout << Color_HighLight <<"(" << RESET;
            cout << exp1.coeff(ep, j);
            cout << Color_HighLight << ")" << RESET;
            if(j!=0 || i!=0) cout << "*" << Color_HighLight << pow(ep,j)*pow(eps,i) << RESET;
            if(j<expr.degree(ep)) cout << " + ";
        }
    }
    if(endlQ) cout << endl;
}

ex Factor(const ex expr_in) {
    ex expr = collect_common_factors(expr_in);
    if(is_a<mul>(expr)) {
        ex ret = 1;
        for(auto item : expr) {
            if(item.has(x(w)) || (item.has(y(w))) || (item.has(z(w)))) ret *= Factor(item);
            else ret *= item;
        }
        return ret;
    } else if(is_a<power>(expr) || expr.match(pow(w1,w2))) {
        return pow(Factor(expr.op(0)), expr.op(1));
    }

    exset xyset;
    expr.find(x(w), xyset);
    expr.find(y(w), xyset);
    expr.find(z(w), xyset);
    expr.find(PL(w), xyset);
    lst xy2s, s2xy;
    for(auto xyi : xyset) {
        symbol txy;
        xy2s.append(xyi==txy);
        s2xy.append(txy==xyi);
    }
    ex expr2 = xyz_pow_simplify(expr);
    expr2 = collect_common_factors(expr2);
    expr2 = expr2.subs(xy2s);
    expr2 = factor(expr2);
    expr2 = expr2.subs(s2xy);
    expr2 = xyz_pow_simplify(expr2);
    expr2 = collect_common_factors(expr2);
    return expr2;
}

ex FactorOutX(const ex expr) {
    exset xset;
    expr.find(x(w), xset);
    lst x2s, s2x;
    for(auto xi : xset) {
        symbol tx;
        x2s.append(xi==tx);
        s2x.append(tx==xi);
    }
    ex expr2 = xyz_pow_simplify(expr);
    expr2 = collect_common_factors(expr2);
    expr2 = expr2.subs(x2s);
    expr2 = factor(expr2);
    expr2 = expr2.subs(s2x);
    expr2 = xyz_pow_simplify(expr2);
    expr2 = collect_common_factors(expr2);
    if(!is_a<mul>(expr2)) return expr2;
    ex ret = 1;
    for(auto item : expr2) {
        if(!item.match(x(w)) && !item.match(pow(x(w1),w2))) ret *= item;
    }
    return ret;
}

ex exp_simplify(const ex expr_in) {
    auto expr = expr_in;
    exmap sub_exp;
    sub_exp[pow(exp(w1),w2)]=exp(w1*w2);
    sub_exp[sqrt(exp(w1))]=exp(w1/2);
    sub_exp[exp(w1)*exp(w2)*w0]=exp(w1+w2)*w0;
    sub_exp[exp(w1)*exp(w2)]=exp(w1+w2);
    while(true) {
        auto expo = expr.subs(sub_exp);
        if(is_zero(expo-expr)) break;
        expr = expo;
    }
    return expr;
}

ex pow_simplify(const ex expr_in) {
    auto expr = expr_in;
    exmap sub_pow;
    sub_pow[pow(pow(w1,w2),w3)] = pow(w1,w2*w3);
    sub_pow[sqrt(pow(w1,w2))] = pow(w1,w2/2);
    sub_pow[pow(sqrt(w1),w2)] = pow(w1,w2/2);
    sub_pow[pow(w1,w2)*pow(w1,w3)*w0] = pow(w1,w2+w3)*w0;
    sub_pow[pow(w1,w2)*w1*w0] = pow(w1,w2+1)*w0;
    sub_pow[pow(w1,w2)/w1*w0] = pow(w1,w2-1)*w0;
    sub_pow[pow(w1,w2)*sqrt(w1)] = pow(w1,w2+1/ex(2));
    sub_pow[pow(w1,w2)/sqrt(w1)] = pow(w1,w2-1/ex(2));
    while(true) {
        auto expo = expr.subs(sub_pow);
        if(is_zero(expo-expr)) break;
        expr = expo;
    }
    return expr;
}

ex xyz_pow_simplify(const ex expr_in) {
    ex expr = expr_in;
    
    //copied from pow_simplify
    exmap sub_pow;
    sub_pow[pow(pow(w1,w2),w3)] = pow(w1,w2*w3);
    sub_pow[sqrt(pow(w1,w2))] = pow(w1,w2/2);
    sub_pow[pow(sqrt(w1),w2)] = pow(w1,w2/2);
    sub_pow[pow(w1,w2)*pow(w1,w3)*w0] = pow(w1,w2+w3)*w0;
    sub_pow[pow(w1,w2)*w1*w0] = pow(w1,w2+1)*w0;
    sub_pow[pow(w1,w2)/w1*w0] = pow(w1,w2-1)*w0;
    sub_pow[pow(w1,w2)*sqrt(w1)] = pow(w1,w2+1/ex(2));
    sub_pow[pow(w1,w2)/sqrt(w1)] = pow(w1,w2-1/ex(2));
    
    sub_pow[pow(x(w1)*w2,w3)] = pow(x(w1),w3)*pow(w2,w3);
    sub_pow[pow(y(w1)*w2,w3)] = pow(y(w1),w3)*pow(w2,w3);
    sub_pow[pow(z(w1)*w2,w3)] = pow(z(w1),w3)*pow(w2,w3);
    while(true) {
        auto expo = expr.subs(sub_pow);
        if(is_zero(expo-expr)) break;
        expr = expo;
    }
    return expr;
}

ex SecDec::PrefactorFIESTA(int nLoop) {
    return  pow(I*pow(Pi,2-ep)*exp(ex(0)-ep*Euler), ex(0)-ex(nLoop));
}


/*-----------------------------------------------------*/
// Refined F-Term
/*-----------------------------------------------------*/
ex SecDec::XRefined(ex const & in_ft) {
    auto ft = Factor(in_ft);
    while(true) {
        auto ft0 = ft;
        if(ft.match(pow(w1, w2))) {
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
    return ft;
}

lst SecDec::XRefined_lst(ex const & in_ft) {
    auto ft = XRefined(in_ft);
    lst ret;
    if(is_a<mul>(ft)) {
        for(auto const &item : ft) ret.append(item);
    } else {
        ret.append(ft);
    }
    return ret;
}

}
